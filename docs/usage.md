# Installation and usage

## Supported inputs

The initial workflow supports local, paired-end Illumina FASTQ files (`.fastq`, `.fq`, optionally `.gz`), an uncompressed FASTA, and optional BED intervals. One row represents one biological sample and one paired library. Multi-lane merging, single-end, BAM/CRAM input, remote object storage input, and cohort genotyping are future features.

The samplesheet header is exactly `sample,fastq_1,fastq_2`. IDs must start with a letter or digit and contain only letters, digits, dots, underscores, and hyphens (maximum 80 characters). FASTQ basenames use the same character set and may also start with an underscore. Paths may contain spaces in parent directories. FASTQ content must agree with its compression suffix. Both mates must have matching read identifiers, record counts, and valid four-line records.

BED uses zero-based, half-open coordinates. Intervals must be nonempty, non-overlapping, within contig bounds, and sorted in FASTA contig order. BED contig names are checked explicitly. These checks cannot establish assembly identity if different assemblies reuse the same contig names: obtain the reference and BED from documented, matching sources.

## Local environment

For Linux x86-64, reproduce the exact package builds used by the Docker image:

```bash
micromamba create -n seqforge-ngs -f environment-linux-64.lock
micromamba activate seqforge-ngs
```

Install Nextflow 24.10.6 separately following the [official installation instructions](https://www.nextflow.io/docs/stable/install.html). Nextflow runs on the host, outside the tool container. Java 17 is used in CI. The readable `environment.yml` specifies top-level versions; the explicit lock pins transitive dependencies as well. Updating the environment requires regenerating the lock and rerunning integration checks.

Native Apple Silicon is not covered by this lock. Use a Linux x86-64 host; emulation has not been tested.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--input` | Required | Samplesheet CSV |
| `--fasta` | Required | Uncompressed reference FASTA |
| `--intervals` | Entire supplied reference | Sorted, merged BED for calling |
| `--outdir` | `results` | Published outputs |
| `--max_cpus` | `8` | Per-process CPU ceiling |
| `--index_memory` | `64 GB` | Reference indexing memory request |
| `--align_memory` | `32 GB` | Alignment memory request |
| `--calling_memory` | `8 GB` | HaplotypeCaller memory request |
| `--min_qual` | `30` | Minimum variant QUAL |
| `--min_depth` | `10` | Minimum sample FORMAT/DP |
| `--container` | `seqforge-ngs:0.1.0-dev` | Locally built image or a supplied image URI |

Human genome indexing can require substantial memory and disk. These resource requests are initial settings, not measured WGS performance guarantees. Sorting uses temporary BAMs in the Nextflow work directory. Whole-genome calling is currently one task per sample and is not scattered. Use a small region for the first human-data evaluation.

## Resume

```bash
nextflow run main.nf -profile test,docker -resume
```

Keep both `work/` and `.nextflow/`; they hold task outputs and cache metadata. Retain the original FASTQs and reference unchanged. Use a new output directory for an unrelated run so stale published files cannot be confused with new results. The report JSON lists only samples participating in its run.

## SLURM and Singularity

The profiles provide a starting configuration; no claim of a tested Redhawk deployment is made. Adapt queue, account, mounts, time limits, and memory to your cluster. Run from the repository so Nextflow can access its modules and scripts.

Build on a machine with Docker, then transfer the image archive to the cluster:

```bash
docker save seqforge-ngs:0.1.0-dev -o seqforge-ngs.tar
```

On the cluster:

```bash
singularity build seqforge-ngs.sif docker-archive://seqforge-ngs.tar
nextflow run main.nf -profile singularity,slurm -c cluster.config \
  --container /shared/containers/seqforge-ngs.sif \
  --input samplesheet.csv --fasta /shared/reference/GRCh38.fa
```

Example `cluster.config` (replace site values):

```groovy
process.queue = 'YOUR_PARTITION'
process.clusterOptions = '--account=YOUR_ACCOUNT'
workDir = '/shared/seqforge/work'
```

For sites using the `apptainer` executable, replace the Singularity engine settings in a local config with `apptainer.enabled = true` and `apptainer.autoMounts = true`, and set `process.container` to your SIF. Apptainer has not yet been exercised in CI.

## Troubleshooting

- **Missing input:** FASTQ paths are relative to the CSV directory; reference/BED CLI paths are relative to the launch directory.
- **Wrong contig name:** `1` and `chr1` are distinct identifiers. Verify the assembly instead of renaming blindly.
- **Insufficient memory:** reduce concurrency or use a larger host; reducing a declared memory request does not reduce the aligner's actual requirements.
- **Task failure:** inspect the reported work directory's `.command.err`, `.command.log`, and `.command.sh`. Correct the cause, then resume.
- **No PASS variants:** check the raw VCF, filter labels, read depth, and mapping. A successfully completed workflow is not proof of adequate data or sensitivity.
