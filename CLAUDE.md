# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Detailed naming, style, shared-code, and output/log-format conventions live in
> `.claude/CLAUDE.md` and are loaded automatically — read them before adding or
> editing a `bin/` script or `modules/*.nf` file. This file covers commands and
> architecture.

## What this is

A containerized Nextflow pipeline (`mg-proc.nf`) for quality-checking,
preprocessing, de novo assembly, and read mapping of shotgun metagenomic
sequencing data (paired-end or single-end). Each pipeline step is a Nextflow
process (`modules/*.nf`) that wraps a standalone script (`bin/*.py` or
`bin/*.R`) and runs in its own per-module Docker image.

## Commands

### Run the pipeline (containerized, via Nextflow)

```bash
# Paired-end, on the bundled test data
nextflow run mg-proc.nf --input_tsv tests/data/samplesheet_pe.tsv --subsample t

# Single-end (reads2 is empty in the samplesheet, so this is inferred automatically)
nextflow run mg-proc.nf --input_tsv tests/data/samplesheet_se.tsv --reformat t --subsample t

# Skip assembly (QC + preprocess only)
nextflow run mg-proc.nf --skip_assembly true

# Full parameter listing
nextflow run mg-proc.nf --help
```

Input is a TSV samplesheet (`--input_tsv`, tab-delimited) with columns
`sample_name`, `reads1`, `reads2`; leave `reads2` empty to mark a sample
single-end (a sheet must be all paired-end or all single-end — mixed sheets
are rejected). `mg-proc.nf` parses it once at startup to build the reads
channel and infer paired-end vs. single-end for the whole run — there is no
separate `--single_end` flag. See `tests/data/samplesheet_pe.tsv` /
`samplesheet_se.tsv` for examples.

Docker images are pulled automatically from `ghcr.io/pereiramemo/mg-proc/*`
(`docker.enabled = true` in `nextflow.config`); no local build is needed unless
a Dockerfile or a `docker/resources/*.requirements.yml` pin changes.

Two params control CPU usage: `--nslots` (threads per tool invocation) and
`--maxForks` (parallel per-sample tasks). Keep `nslots * maxForks` at or below
the host's physical core count — it's the local executor's CPU budget
(`executor.cpus` in `nextflow.config`), and `MODULE_1_2_QUALITY_CHECK` alone
reserves that entire budget since it aggregates all samples at once.

### Run the `bin/` scripts directly (no containers, no Nextflow)

```bash
bash tests/run_tests.sh
```

Script-level smoke test for `2-preprocess.py`, `3-assembly-and-map.py`, and
`1.1-quality-check.py` (SE and PE) plus a `py_compile` syntax check. Requires
fastp, bbmap, seqtk, pear, pigz, megahit, bwa, samtools, picard, and
matplotlib on `PATH` (e.g. a conda/mamba env) — the scripts shell out to these
tools directly. Each script can also be invoked standalone, e.g.:

```bash
bin/2-preprocess.py --reads R1.fastq.gz --reads2 R2.fastq.gz --sample_name s1 --output_dir out --overwrite t
```

There is no dedicated test runner for the R script (`1.2-quality-check.R`);
exercise it through `nextflow run` or by invoking it directly with `Rscript`.

### Syntax-check a single script

```bash
python -m py_compile bin/2-preprocess.py
```

### Build/publish Docker images

Only needed after changing a `docker/*.Dockerfile` or a
`docker/resources/*.requirements.yml` pin — run from the repo root:

```bash
bash docker/dockerbuild_commands.sh                          # build + tag :latest locally
PUSH=1 VERSION=v1.0.0 bash docker/dockerbuild_commands.sh    # build, tag :latest and :v1.0.0, push
```

Requires `docker login ghcr.io` first when pushing. Newly pushed packages are
private by default and must be made public on GitHub for anonymous pulls.

## Architecture

### Module ↔ script ↔ container correspondence

Each pipeline step is three files that must stay in lockstep, tied together by
the shared base name (see `.claude/CLAUDE.md` for the exact naming scheme):

- `bin/<name>.py` or `.R` — the actual logic; runs standalone, invoked by name
  because Nextflow stages the whole `bin/` directory onto `PATH`.
- `modules/<name>.nf` — a Nextflow `process` wrapping that script: declares
  the container image, `publishDir`, inputs/outputs, and builds the CLI
  invocation from `params.*` / `task.cpus`.
- `docker/<name>.Dockerfile` + `docker/resources/<name>.requirements.yml` —
  the per-module image (micromamba-based) with pinned conda/mamba deps.

`mg-proc.nf` is the workflow entry point: it builds the reads channel (PE via
`channel.fromFilePairs`, SE via `channel.fromPath`), then wires the four
modules together. `nextflow.config` holds every `params.*` default and the
executor CPU budget.

### Pipeline flow

```
reads ─┬─> MODULE_1_1_QUALITY_CHECK      (per-sample fastp report, diagnostic, always runs)
       ├─> MODULE_1_2_QUALITY_CHECK      (aggregate R/ShortRead/DADA2 plots, diagnostic, always runs)
       └─> MODULE_2_PREPROCESS ──(QC-trimmed reads: *_qc-02.fastq[.gz])──> MODULE_3_ASSEMBLY_AND_MAP
```

`MODULE_2_PREPROCESS` emits the reads that feed assembly via its `qc_reads`
output channel — paired-end `*_R{1,2}_qc-02.fastq`, single-end
`*_se_qc-02.fastq`. `MODULE_3_ASSEMBLY_AND_MAP` is skipped with
`--skip_assembly true`. Single-end vs. paired-end is a pipeline-wide switch inferred from
`--input_tsv` (empty `reads2` column => single-end); `mg-proc.nf` computes it
once and passes it as an explicit `single_end` process input to every module,
which branches on it internally rather than having separate SE/PE modules.

### Shared helpers (`bin/utils.py` / `bin/utils.R`)

Every script imports/sources its language's helper module for: the
`log`/`log_warn`/`log_error` console+buffer helpers, `derive_sample_name`, and
`build_log` (assembles the standardized general-info log header, see
`.claude/CLAUDE.md`'s "Log format" section). `utils.py` additionally provides
FASTQ/FASTA counting helpers (`count_fastq`, `count_fasta`, `mean_length`)
used for the `stats/` TSVs, and both modules provide `decompress_or_link(src,
dst)` — normalizes a read file to plain text (decompress gzip/bzip2, or
symlink when already plain), needed before handing reads to a tool that
doesn't understand bzip2 (e.g. fastp) or that operates on a whole directory
(e.g. `1.2-quality-check.R`'s ShortRead calls, which stage a decompressed copy
of every matched file before running `qa()`/`FastqSampler()`). The two modules
are kept in sync by hand so Python and R scripts emit identical log layouts —
when changing one, check whether the other needs the same change.

### Script internal structure

Each `bin/*.py` script follows the same shape: `parse_args()` → validate
inputs/required tools (`shutil.which`) → create `output/`, `logs/`, `stats/`
under `--output_dir` → run a sequence of numbered steps (each tool invocation
goes through a local `run_tool`/`run_redirect` wrapper that also appends to an
accumulated `tool_log` list) → write `stats/*.tsv` → write `logs/*.log` via
`build_log(...)`. Intermediate FASTQ filenames encode the pipeline stage as a
trailing `-NN` suffix (e.g. `_R1_at-01.fastq`, `_assembled_qc-03.fastq`); this
ordering is parsed back out (`parse_order`) when building `stats.tsv`.

### Test data vs. real runs

`tests/data/` holds small bundled FASTQ sets (`SRR12479690` paired-end,
`SRR4831661` single-end) used as pipeline defaults and by `tests/run_tests.sh`,
plus the `samplesheet_pe.tsv` / `samplesheet_se.tsv` samplesheets that point at
them. `tests/data_samo/` and `tests/run_tests_samo.sh` are a separate, ad hoc
local dataset/run script (hostname-specific absolute paths in
`tests/data_samo/samplesheet.tsv`) — not part of the standard test suite.
