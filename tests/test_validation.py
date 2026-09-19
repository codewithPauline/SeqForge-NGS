import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from validate_fastq import validate_pair
from validate_reference import bed_intervals, fasta_lengths
from make_report import render_report, samtools_stats, variant_stats


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, content):
        path = self.root / name
        path.write_text(content)
        return path

    def pair(self, one='@read/1\nACGT\n+\nIIII\n', two='@read/2\nTGCA\n+\nIIII\n'):
        return self.write('r1.fastq', one), self.write('r2.fastq', two)

    def test_valid_mates(self):
        self.assertEqual(validate_pair(*self.pair()), 1)

    def test_gzip_mates(self):
        paths = []
        for mate in (1, 2):
            path = self.root / f'r{mate}.fastq.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write(f'@read {mate}:N:0:ACGT\nACGT\n+\nIIII\n')
            paths.append(path)
        self.assertEqual(validate_pair(*paths), 1)

    def test_unequal_record_counts(self):
        one, two = self.pair()
        two.write_text(two.read_text() * 2)
        with self.assertRaisesRegex(ValueError, 'Unequal'):
            validate_pair(one, two)

    def test_mismatched_ids(self):
        with self.assertRaisesRegex(ValueError, 'Mismatched'):
            validate_pair(*self.pair(two='@other/2\nACGT\n+\nIIII\n'))

    def test_truncated_quality(self):
        with self.assertRaisesRegex(ValueError, 'length mismatch'):
            validate_pair(*self.pair(two='@read/2\nACGT\n+\nIII\n'))

    def test_wrong_orientation(self):
        with self.assertRaisesRegex(ValueError, 'orientation'):
            validate_pair(*self.pair(two='@read/1\nACGT\n+\nIIII\n'))

    def test_empty_fastq(self):
        with self.assertRaisesRegex(ValueError, 'no reads'):
            validate_pair(*self.pair('', ''))

    def test_invalid_bases(self):
        with self.assertRaisesRegex(ValueError, 'DNA'):
            validate_pair(*self.pair(one='@read/1\nACGX\n+\nIIII\n'))

    def test_invalid_quality(self):
        with self.assertRaisesRegex(ValueError, 'quality encoding'):
            validate_pair(*self.pair(one='@read/1\nACGT\n+\nIII \n'))

    def test_mislabeled_gzip(self):
        one, two = self.pair()
        one = one.rename(self.root / 'r1.fastq.gz')
        with self.assertRaisesRegex(ValueError, 'gzip content'):
            validate_pair(one, two)

    def test_reference_multiline_and_case(self):
        reference = self.write('ref.fa', '>chr1 comment\nACGT\nacgt\n>chr2\nNNNN\n')
        self.assertEqual(fasta_lengths(reference), {'chr1': 8, 'chr2': 4})

    def test_duplicate_contig(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            fasta_lengths(self.write('ref.fa', '>chr1\nACGT\n>chr1\nTGCA\n'))

    def test_empty_contig(self):
        with self.assertRaisesRegex(ValueError, 'empty contig'):
            fasta_lengths(self.write('ref.fa', '>chr1\n>chr2\nACGT\n'))

    def test_sequence_before_header(self):
        with self.assertRaisesRegex(ValueError, 'Invalid FASTA sequence'):
            fasta_lengths(self.write('ref.fa', 'ACGT\n>chr1\nACGT\n'))

    def test_valid_bed_zero_based_half_open(self):
        bed = self.write('target.bed', 'chr1\t0\t4\nchr1\t4\t8\n')
        self.assertEqual(bed_intervals(bed, {'chr1': 8}), [('chr1', 0, 4), ('chr1', 4, 8)])

    def test_reference_mismatch(self):
        with self.assertRaisesRegex(ValueError, 'not found'):
            bed_intervals(self.write('target.bed', '1\t0\t4\n'), {'chr1': 8})

    def test_out_of_bounds(self):
        for start, end in [(-1, 4), (1, 9), (4, 4)]:
            with self.subTest(start=start, end=end), self.assertRaisesRegex(ValueError, 'bounds'):
                bed_intervals(self.write('target.bed', f'chr1\t{start}\t{end}\n'), {'chr1': 8})

    def test_overlapping_or_unsorted_bed(self):
        for text in ['chr1\t0\t5\nchr1\t4\t8\n', 'chr2\t0\t3\nchr1\t0\t3\n']:
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, 'sorted'):
                bed_intervals(self.write('target.bed', text), {'chr1': 8, 'chr2': 8})

    def test_missing_metrics_fail_loudly(self):
        for parser in (samtools_stats, variant_stats):
            with self.subTest(parser=parser), self.assertRaisesRegex(ValueError, 'Missing required'):
                parser(self.write('empty.stats', '# no metrics\n'))

    def test_report_escapes_metadata(self):
        report = render_report({'samples': [], 'version': '<script>alert(1)</script>'})
        self.assertNotIn('<script>', report)
        self.assertIn('&lt;script&gt;', report)


if __name__ == '__main__':
    unittest.main()
