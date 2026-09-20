import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from benchmark_io import download, intersect_bed, inventory, verify, REQUIRED
from benchmark_report import read_roc


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def roc(self, baseline, calls, rows):
        path = self.root / 'test.tsv.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write(f'#total baseline variants: {baseline}\n#total call variants: {calls}\n')
            handle.write('#score\ttrue_positives_baseline\tfalse_positives\ttrue_positives_call\tfalse_negatives\tprecision\tsensitivity\tf_measure\n')
            handle.write(rows)
        return path

    def test_reports_all_calls_instead_of_best_f1(self):
        result = read_roc(self.roc(2, 2, '100\t1\t0\t1\t1\t1.0\t0.5\t0.6667\n30\t1\t1\t1\t1\t0.5\t0.5\t0.5\n'))
        self.assertEqual(result['precision'], 0.5)
        self.assertEqual(result['f1'], 0.5)
        self.assertEqual(result['lowest_included_qual'], '30')

    def test_no_calls_does_not_invent_perfect_precision(self):
        result = read_roc(self.roc(2, 0, ''))
        self.assertIsNone(result['precision'])
        self.assertIsNone(result['f1'])
        self.assertEqual(result['recall'], 0)

    def test_empty_truth_and_query(self):
        result = read_roc(self.roc(0, 0, ''))
        self.assertIsNone(result['recall'])

    def test_incomplete_roc_fails(self):
        with self.assertRaisesRegex(ValueError, 'account'):
            read_roc(self.roc(2, 2, '100\t1\t0\t1\t1\t1\t0.5\t0.6667\n'))

    def test_rescaled_truth_and_query_counts_remain_distinct(self):
        result = read_roc(self.roc(3, 2, '40\t3\t0\t2\t0\t1\t1\t1\n'))
        self.assertEqual((result['tp_baseline'], result['tp_call']), (3, 2))

    def test_interval_clipping_merging_and_contig_selection(self):
        path = self.root / 'regions.bed'
        path.write_text('chr1\t0\t20\nchr1\t18\t25\nchr1\t25\t40\nchr2\t0\t100\n')
        self.assertEqual(intersect_bed(path, 'chr1', 10, 30), [(10, 30)])

    def test_empty_confident_intersection_fails(self):
        path = self.root / 'regions.bed'
        path.write_text('chr1\t0\t10\n')
        with self.assertRaisesRegex(ValueError, 'no overlap'):
            intersect_bed(path, 'chr1', 10, 20)

    def asset(self, data=b'abc'):
        return {'filename': 'test.bin', 'url': 'https://example.org/test.bin', 'bytes': len(data), 'md5': hashlib.md5(data).hexdigest()}

    def test_download_checks_provider_checksum(self):
        with patch('urllib.request.urlopen', return_value=io.BytesIO(b'abc')):
            self.assertEqual(download(self.asset(), self.root).read_bytes(), b'abc')

    def test_download_rejects_wrong_bytes(self):
        with patch('urllib.request.urlopen', return_value=io.BytesIO(b'xyz')), self.assertRaisesRegex(ValueError, 'MD5'):
            download(self.asset(), self.root)
        self.assertFalse((self.root / 'test.bin').exists())
        self.assertFalse((self.root / 'test.bin.part').exists())

    def test_download_rejects_oversize(self):
        with patch('urllib.request.urlopen', return_value=io.BytesIO(b'abcd')), self.assertRaisesRegex(ValueError, 'exceeds'):
            download(self.asset(), self.root)

    def test_asset_path_traversal_rejected(self):
        asset = self.asset()
        asset['filename'] = '../outside.bin'
        with self.assertRaisesRegex(ValueError, 'basename'):
            download(asset, self.root)

    def test_prepared_file_mutation_detected(self):
        for name in REQUIRED:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('abc')
        manifest = {'schema_version': 1, 'sample': 'TEST', 'padding_bases': 10, 'min_qual': 30, 'min_depth': 10,
                    'scoring_region': {'chrom': 'chr1', 'start': 1, 'end': 100, 'coordinate_system': '0-based half-open'}}
        (self.root / 'prepared.json').write_text(json.dumps({'schema_version': 1, 'manifest': manifest, 'files': inventory(self.root)}))
        self.assertTrue(verify(self.root)['verified'])
        (self.root / 'reads.bam').write_text('xyz')
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            verify(self.root)


if __name__ == '__main__':
    unittest.main()
