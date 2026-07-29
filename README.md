# Metagenomic Processing Pipeline

This repository provides a containerized [Nextflow](https://www.nextflow.io/)
pipeline for quality checking, preprocessing, de novo assembly, and read mapping
of shotgun metagenomic sequencing data (paired-end or single-end). Each step is a
Nextflow process (under `modules/`) that wraps a Python or R script in `bin/` and
runs in its own Docker image.

The workflow runs end to end: quality-check diagnostics, then preprocessing, then
assembly + mapping of the preprocessed reads. Assembly can be skipped with
`--skip_assembly`.

## Repository structure

```
.
├── LICENSE                                 # License file
├── README.md                               # This file
├── mg-proc.nf                              # Nextflow workflow entry point
├── nextflow.config                         # Nextflow parameters and Docker settings
├── bin/                                    # Step scripts (auto-staged onto PATH)
│   ├── 1.1-quality-check.py                # fastp QC report (per sample)
│   ├── 1.2-quality-check.R                 # Comparative QC plots (all samples)
│   ├── 2-preprocess.py                     # Preprocessing pipeline (per sample)
│   ├── 3-assembly-and-map.py               # Assembly + mapping (per sample)
│   ├── utils.py                            # Shared Python helpers
│   └── utils.R                             # Shared R helpers
├── modules/                                # Nextflow process definitions (*.nf)
├── docker/                                 # Per-module Dockerfiles + build script
│   ├── *.Dockerfile
│   ├── dockerbuild_commands.sh
│   └── resources/*.requirements.yml
└── tests/                                  # Bundled test data + smoke test
```

## Installation

The pipeline runs entirely in containers, so the only prerequisites are:

- **Java** 11 or later (required by Nextflow)
- **[Nextflow](https://www.nextflow.io/)** 23.04 or later:
  ```bash
  curl -s https://get.nextflow.io | bash
  sudo mv nextflow /usr/local/bin/    # or any directory on your PATH
  ```
- **[Docker](https://docs.docker.com/get-docker/)** (the daemon must be running; the
  invoking user must be able to run `docker`)

Then clone the repository:

```bash
git clone https://github.com/pereiramemo/MG-Proc.git
cd MG-Proc
```

The per-module Docker images are published at `ghcr.io/pereiramemo/mg-proc/*`, and
Nextflow pulls them automatically on the first `nextflow run` (`docker.enabled = true`
in `nextflow.config`). You only need to build images yourself if you change a
Dockerfile or a pinned dependency — see
[Building & publishing the images](#building--publishing-the-images).

## Pipeline steps

| Module | Script | Purpose |
|--------|--------|---------|
| `MODULE_1_1_QUALITY_CHECK`   | `1.1-quality-check.py`   | fastp QC report, report-only (diagnostic, always runs) |
| `MODULE_1_2_QUALITY_CHECK`   | `1.2-quality-check.R`    | Comparative QC plots: mean-quality vs read count, count histograms, PhiX % (diagnostic, always runs) |
| `MODULE_2_PREPROCESS`        | `2-preprocess.py`        | Reformat/repair, subsample, adapter trim, quality trim, PE merge, FASTA conversion |
| `MODULE_3_ASSEMBLY_AND_MAP`  | `3-assembly-and-map.py`  | MEGAHIT assembly + BWA-MEM mapping + optional Picard dedup |

Each step writes a standardized layout under its publish directory: `output/`
(main results), `logs/` (a log file with a general-info header followed by any
third-party tool output), and `stats/` (TSV statistics). The output, logging, and
naming conventions are documented in `.claude/CLAUDE.md`.

## Workflow

```
reads ─┬─> MODULE_1_1_QUALITY_CHECK      (per sample, diagnostic)
       ├─> MODULE_1_2_QUALITY_CHECK      (all samples, diagnostic)
       └─> MODULE_2_PREPROCESS ──(QC-trimmed reads)──> MODULE_3_ASSEMBLY_AND_MAP
```

`MODULE_1_1_QUALITY_CHECK` (fastp) and `MODULE_1_2_QUALITY_CHECK` (comparative plots)
are diagnostics that always run on the raw reads. `MODULE_2_PREPROCESS` runs per
sample; the quality-trimmed reads it emits (paired-end `*_qc-02.fastq`, single-end
`*_se_qc-02.fastq`) feed `MODULE_3_ASSEMBLY_AND_MAP`, which is skipped when
`--skip_assembly true`.

Paired-end (default) and single-end reads are both supported: with `--single_end true`
the reads channel is built with `channel.fromPath` (using `--se_reads_pattern`) instead
of `channel.fromFilePairs`, and every module runs in single-end mode.

## Run

```bash
# Paired-end run on the bundled test data (reads_pattern matches the SRA _1/_2 names)
nextflow run mg-proc.nf \
    --input_dir tests/data/SRR12479690 \
    --reads_pattern '*_{1,2}.fastq.gz' \
    --subsample t

# Single-end run
nextflow run mg-proc.nf \
    --single_end true \
    --input_dir tests/data/SRR4831661 \
    --se_reads_pattern '*.fastq' \
    --reformat t --subsample t

# QC + preprocess only (no assembly)
nextflow run mg-proc.nf --skip_assembly true

# On your own data
nextflow run mg-proc.nf \
    --input_dir     /path/to/fastq \
    --reads_pattern '*_R{1,2}_001.fastq.gz' \
    --output_dir    /path/to/results \
    --trim_adapters t \
    --nslots        16

# Full parameter listing
nextflow run mg-proc.nf --help
```

## Parameters

All parameters have defaults in `nextflow.config` and can be overridden on the command
line (e.g. `--nslots 16`). The full list (output of `nextflow run mg-proc.nf --help`):

```text
General:
  --input_dir         DIR   Input directory with FASTQ files (default: ./tests/data)
  --reads_pattern     STR   Paired-end glob for fromFilePairs (default: *_R{1,2}_001.fastq.gz)
  --se_reads_pattern  STR   Single-end glob when --single_end true (default: *.fastq.gz)
  --single_end        BOOL  Process single-end reads (default: false)
  --output_dir        DIR   Output directory (default: ./output_nf)
  --nslots            INT   CPU threads per tool (default: 12)
  --maxForks          INT   Max parallel process instances (default: 3)
  --full_output       BOOL  Publish all module outputs (default: true)
  --skip_assembly     BOOL  Skip MODULE_3_ASSEMBLY_AND_MAP (default: false)
  --container_tag     STR   Tag of the ghcr.io/pereiramemo/mg-proc/* images (default: latest)

MODULE_1_1_QUALITY_CHECK — fastp QC report:
  --qc_min_length             INT  Minimum read length, reporting only (default: 50)
  --qualified_quality_phred   INT  Qualified base quality, reporting only (default: 20)
  --unqualified_percent_limit INT  Max unqualified base percent (default: 40)
  --disable_adapter_trimming  STR  Disable adapter trimming in report, t/f (default: t)

MODULE_1_2_QUALITY_CHECK — comparative QC plots:
  (finds files via --reads_pattern / --se_reads_pattern from General)
  --qc_sample_size  INT   Reads subsampled per file for QC (default: 10000)

MODULE_2_PREPROCESS — preprocessing:
  --reformat        STR  Reformat FASTQ with reformat.sh, t/f (default: f)
  --repair          STR  Reformat + repair FASTQ (PE), t/f (default: f)
  --subsample       STR  Subsample to 10k reads, t/f (default: f)
  --trim_adapters   STR  Remove adapters with BBDuk, t/f (default: f)
  --output_merged   STR  Merge PE reads, t/f (default: t)
  --merger          STR  pear | bbmerge (default: pear)
  --min_overlap     INT  Minimum PE overlap for PEAR (default: 10)
  --pvalue          NUM  p-value for PEAR (default: 0.01)
  --min_length      INT  Minimum read length after trimming (default: 75)
  --min_qual        INT  Quality trim threshold (default: 20)
  --seed            INT  Random seed for subsampling (default: 123)
  --clean           STR  Remove intermediates, t/f (default: f)
  --compress        STR  Compress outputs with pigz, t/f (default: f)

MODULE_3_ASSEMBLY_AND_MAP — assembly + mapping:
  --assem_preset      STR  MEGAHIT preset (default: meta-sensitive)
  --min_contig_length INT  Minimum contig length to keep (default: 250)
  --remove_duplicates STR  Remove PCR duplicates with Picard, t/f (default: f)
  --contigs           STR  Pre-assembled contigs FASTA (default: none)
  --assem_dir         STR  Directory with previous assemblies (default: none)
```

## Threads and parallelism

Two parameters control CPU usage:

- `--nslots` — threads given to one tool invocation (a single task).
- `--maxForks` — how many per-sample tasks run at the same time.

Per-sample modules (`MODULE_1_1_QUALITY_CHECK`, `MODULE_2_PREPROCESS`,
`MODULE_3_ASSEMBLY_AND_MAP`) reserve `cpus = nslots` each, so up to `maxForks` run at
once. The aggregate `MODULE_1_2_QUALITY_CHECK` reads every sample together (its input is
`.collect()`ed in `mg-proc.nf`) and reserves the whole budget `nslots * maxForks`, so it
runs alone with all threads. This is enforced by the local-executor budget
`executor.cpus = nslots * maxForks` in `nextflow.config`.

> **Keep `nslots * maxForks` at or below the machine's physical core count.** The executor
> budget is pinned to that product, so a larger value oversubscribes the CPUs — Nextflow
> will not clamp it for you. Example on a 48-core host:
> `nextflow run mg-proc.nf --nslots 16 --maxForks 3` (16 × 3 = 48).

## Building & publishing the images

End users do **not** need this section — the published images pull automatically. It is
only for rebuilding and republishing after changing a Dockerfile or a pinned dependency
in `docker/resources/*.requirements.yml`. The images are built from the per-module
Dockerfiles in `docker/` by `docker/dockerbuild_commands.sh` (run from the repository
root):

```bash
# Build + tag :latest locally
bash docker/dockerbuild_commands.sh

# Build, tag with a version, and push (:latest and :v1.0.0) to the registry
echo "$GHCR_PAT" | docker login ghcr.io -u pereiramemo --password-stdin   # PAT needs write:packages
PUSH=1 VERSION=v1.0.0 bash docker/dockerbuild_commands.sh
```

The script honours two environment variables: `VERSION` (adds an extra immutable tag
alongside `:latest`) and `PUSH=1` (pushes after building). Newly pushed packages are
**private by default**; make each one public (GitHub → **Packages** → **Package settings**
→ **Change visibility** → **Public**) so machines can pull them anonymously.

### Reproducible installs (image version pinning)

Every module pulls `ghcr.io/pereiramemo/mg-proc/<module>:${params.container_tag}`.
`container_tag` defaults to `latest`; for a reproducible install, pin a published version,
either per run (`nextflow run mg-proc.nf --container_tag v1.0.0`) or by changing the default
in `nextflow.config`. A pinned tag must already be published, or the pull fails.

## Dependencies

Dependencies are pinned per module in `docker/resources/*.requirements.yml` and built
into the per-module images — there is nothing to install manually beyond Nextflow and
Docker.

| Tool | Purpose |
|------|---------|
| [fastp](https://github.com/OpenGene/fastp) | Read quality control (QC report) |
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) | reformat.sh, repair.sh, BBDuk, BBMerge (reformat/repair, adapter & quality trim, merging) |
| [seqtk](https://github.com/lh3/seqtk) | Subsampling and FASTQ→FASTA conversion |
| [PEAR](https://cme.h-its.org/exelixis/web/software/pear) | Paired-end read merging |
| [pigz](https://zlib.net/pigz/) | Parallel gzip compression |
| [MEGAHIT](https://github.com/voutcn/megahit) | De novo metagenomic assembly |
| [BWA](https://github.com/lh3/bwa) | Read mapping |
| [SAMtools](http://www.htslib.org/) | SAM/BAM manipulation |
| [Picard](https://broadinstitute.github.io/picard/) | PCR duplicate removal |
| [R](https://www.r-project.org/) + [tidyverse](https://www.tidyverse.org/) / [ShortRead](https://bioconductor.org/packages/ShortRead/) / [DADA2](https://benjjneb.github.io/dada2/) | Comparative QC plots and PhiX detection |
| [Python 3](https://www.python.org/) | Step scripts and orchestration helpers |

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
