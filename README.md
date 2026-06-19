# Metagenomic pipelines
This repository contains scripts for quality checking, preprocessing, assembling, and mapping metagenomic data.

# Repository structure

```
.
├── LICENSE
├── README.md
├── environment.yml                          # Conda environment specification
└── modules/                                 # Pipeline modules
    ├── 1.1-quality_check_fastp.sh           # Quality check using fastp
    ├── 1.2-quality_check.R                  # Quality check plots using R
    ├── 2-preprocess_pipeline.sh             # Main preprocessing pipeline script
    ├── 3-assembly_and_map_pipeline.sh       # De novo assembly and read mapping
    ├── conf.sh                              # Configuration file with tool paths
    └── resources/                           # Additional resources and scripts
        ├── fq2fa.sh                         # FASTQ to FASTA conversion script
        └── plots.R                          # R script for generating statistics plots
```

# **Installation instructions**

**Clone the repository**:  
```bash
git clone https://github.com/pereiramemo/MG-Proc.git
cd MG-Proc
```

**Installation with mamba**:  
All dependencies can be installed using mamba (or conda). 

First, check if mamba is installed:
```bash
command -v mamba
```

If mamba is not installed, you can install it via:
- **Miniforge** (recommended): [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)
- **Mambaforge**: [https://github.com/conda-forge/miniforge#mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
- **Or install mamba into an existing conda installation**:
  ```bash
  conda install -n base -c conda-forge mamba
  ```

Once mamba is installed, create a new environment using the provided `environment.yml` file:
```bash
mamba env create -f environment.yml
```

Then activate the environment:
```bash
mamba activate MG-Proc
```

# **How to use**

## 1.1-quality_check_fastp.sh

[1.1-quality_check_fastp.sh](modules/1.1-quality_check_fastp.sh): Quality assessment of raw Illumina reads (paired-end or single-end) using fastp. Runs in report-only mode — input files are never modified.

### Main analysis
- Batch quality assessment across all samples in an input directory
- Paired-end and single-end read support
- Per-sample HTML and JSON reports with base quality, GC content, duplication, and length distribution
- Aggregated summary statistics table across all samples

### Output files
- `reports/SAMPLE_fastp.html`: Interactive HTML quality report per sample
- `reports/SAMPLE_fastp.json`: Machine-readable JSON quality report per sample
- `reports/SAMPLE_fastp.log`: fastp run log per sample
- `stats/summary.tsv`: Tab-separated table with read counts and Q20/Q30 base percentages for all samples
- `summary_report.txt`: Plain-text run summary including parameters and aggregated statistics

### Help
```
Usage: 1.1-quality_check_fastp.sh <options>

Options:
    --help
        Print this help message and exit

    --input_dir=CHAR
        Directory containing input FASTQ files (required)

    --output_dir=CHAR
        Directory to output generated data (QC reports and plots) (required)

    --single_end=t|f
        Process single-end reads instead of paired-end [default=f]

    --r1_pattern=CHAR
        Pattern for R1 FASTQ files, or single-end files when --single_end=t
        [default=_R1_001.fastq.gz]

    --r2_pattern=CHAR
        Pattern for R2 FASTQ files (ignored when --single_end=t)
        [default=_R2_001.fastq.gz]

    --nslots=NUM
        Number of threads to use [default=12]

    --min_length=NUM
        Minimum read length filter (only for reporting) [default=50]

    --qualified_quality_phred=NUM
        Minimum quality value for qualified base (Phred score, only for reporting) [default=20]

    --unqualified_percent_limit=NUM
        Maximum percent of unqualified bases allowed (only for reporting) [default=40]

    --disable_adapter_trimming=t|f
        Disable adapter trimming in report [default=t]

    --html_report=t|f
        Generate HTML report [default=t]

    --json_report=t|f
        Generate JSON report [default=t]

    --overwrite=t|f
        Overwrite previous directory [default=f]
```

### Example usage
```bash
./modules/1.1-quality_check_fastp.sh \
    --input_dir=raw_data/ \
    --output_dir=results/qc_reports
```

---

## 1.2-quality_check.R

[1.2-quality_check.R](modules/1.2-quality_check.R): R script that generates summary quality plots from raw Illumina reads (paired-end or single-end). Reads FASTQ files directly and produces publication-ready PNG figures.

### Main analysis
- Calculation of mean quality scores per sample for R1 (and R2 in paired-end mode)
- Scatter plot of mean quality score versus read count
- Histogram of read count distributions (linear and log scale)
- Detection and quantification of PhiX contamination per sample

### Output files
- `r1_mean_q_vs_nseq.png`: Scatter plot of mean R1 quality score vs. number of reads per sample
- `r2_mean_q_vs_nseq.png`: Scatter plot of mean R2 quality score vs. number of reads per sample (paired-end only)
- `samples_hist.png`: Histogram of read counts across samples
- `samples_hist_log.png`: Log-scale histogram of read counts across samples
- `samples_perc_phix_barplot.png`: Bar plot of estimated PhiX contamination percentage per sample

### Help
```
Usage: 1.2-quality_check.R [options]

Options:
        --input_dir=CHARACTER
                Input directory with FASTQ files

        --output_dir=CHARACTER
                Output directory for plots

        --nslots=INTEGER
                Number of threads to use [default=12]

        --single_end=LOGICAL
                Process single-end reads instead of paired-end [default=FALSE]

        --r1_pattern=CHARACTER
                Pattern for R1 FASTQ files, or single-end files when
                --single_end=TRUE [default=R1_001.fastq.gz]

        --r2_pattern=CHARACTER
                Pattern for R2 FASTQ files (ignored when --single_end=TRUE)
                [default=R2_001.fastq.gz]

        --overwrite=LOGICAL
                Overwrite previous output [default=FALSE]

        -h, --help
                Show this help message and exit
```

### Example usage
```bash
Rscript modules/1.2-quality_check.R \
    --input_dir=raw_data/ \
    --output_dir=results/qc_plots
```

---

## 2-preprocess_pipeline.sh

[2-preprocess_pipeline.sh](modules/2-preprocess_pipeline.sh): Preprocessing pipeline for raw Illumina metagenomic reads (paired-end or single-end). Produces a quality-trimmed FASTA file ready for downstream analyses such as assembly and mapping.

### Main analysis
- Optional FASTQ reformat (`--reformat=t`): runs `reformat.sh tossbrokenreads=t` to fix malformed records and SRA-format extended `+` headers (SE and PE)
- Optional FASTQ repair (`--repair=t`): runs `reformat.sh` first, then `repair.sh` to restore mate pairing (PE only; SE gets reformat only)
- Optional adapter trimming with BBDuk
- Optional subsampling to 10,000 reads for rapid testing
- Paired-end read merging with PEAR or BBMerge (PE only)
- Quality trimming of merged reads, unmerged PE reads, or single-end reads with BBDuk
- FASTQ-to-FASTA conversion of quality-trimmed output
- Read count and mean length statistics at each processing step, with optional plots

### Output files
- `SAMPLE_workable.fasta` (or `.fasta.gz`): Quality-trimmed FASTA file ready for downstream analyses — merged reads (PE) or quality-trimmed reads (SE)
- `stats.tsv`: Tab-separated table with read counts and mean lengths at each pipeline step; includes one row per intermediate FASTQ file, including singletons discarded by `repair.sh` when `--repair=t`
- `stats_plots.png`: Plot of read counts and mean lengths across processing steps (when `--plot=t`)
- `SAMPLE_R1_qc-02.fastq` / `SAMPLE_R2_qc-02.fastq`: Quality-trimmed paired-end FASTQ files (when `--output_pe=t`)
- `SAMPLE_assembled_qc-03.fasta`: Quality-trimmed merged reads in FASTA (PE, when `--output_merged=t`)
- `SAMPLE_unassembled_R1_qc-03.fasta` / `SAMPLE_unassembled_R2_qc-03.fasta`: Quality-trimmed unmerged reads in FASTA (PE, when `--output_merged=t`)
- `SAMPLE_singletons_repaired-00.fastq`: PE reads that lost their mate during `repair.sh` and could not be re-paired (PE only, when `--repair=t`); included in `stats.tsv` as a record of discarded reads

### Help
```
Usage: 2-preprocess_pipeline.sh [OPTIONS]

Required:
  --reads FILE          R1 file (or single-end file when --single_end=t)
  --output_dir DIR      Output directory

Paired-end only (ignored when --single_end=t):
  --reads2 FILE         R2 file
  --merger STR          pear|bbmerge (default pear)
  --min_overlap NUM     Minimum PE overlap for PEAR (default 10)
  --output_pe t|f       Output QC'ed paired-end reads (default f)
  --output_merged t|f   Output merged QC'ed reads (default t)
  --pvalue NUM          p-value for PEAR (default 0.01)

Optional:
  --single_end t|f      Process as single-end reads (default f)
  --reformat t|f        Reformat FASTQ with reformat.sh (default f)
                        Fixes malformed records and SRA extended + headers
  --repair t|f          Reformat + repair FASTQ (default f)
                        Runs reformat.sh then repair.sh (repair.sh PE only)
  --clean t|f           Remove intermediates (default f)
  --compress t|f        Compress outputs with pigz (default f)
  --min_length NUM      Minimum read length after trimming (default 75)
  --min_qual NUM        Quality trim threshold (default 20)
  --nslots NUM          Threads (default 12)
  --overwrite t|f       Replace existing output dir (default f)
  --plot t|f            Produce QC plots (default f)
  --sample_name STR     Name prefix (default metagenomex)
  --seed NUM            Random seed for subsampling (default 123)
  --subsample t|f       Subsample to 10k reads (default f)
  --trim_adapters t|f   Remove adapters (default f)
  --help                Show this help
```

### Example usage
```bash
# Paired-end with adapter trimming and merging
./modules/2-preprocess_pipeline.sh \
    --reads sample_R1.fastq.gz \
    --reads2 sample_R2.fastq.gz \
    --output_dir results/sample1_preproc \
    --trim_adapters=t \
    --output_merged=t \
    --sample_name=sample1

# Single-end with repair (e.g. SRA data)
./modules/2-preprocess_pipeline.sh \
    --reads sample_SE.fastq \
    --single_end=t \
    --repair=t \
    --output_dir results/sample1_preproc \
    --sample_name=sample1
```

---

## 3-assembly_and_map_pipeline.sh

[3-assembly_and_map_pipeline.sh](modules/3-assembly_and_map_pipeline.sh): De novo assembly and read mapping pipeline for metagenomic samples. Supports paired-end and single-end reads, and can use pre-assembled contigs instead of running assembly.

### Main analysis
- De novo assembly with MEGAHIT using paired-end (`-1`/`-2`) or single-end (`-r`) mode; alternatively, accepts pre-assembled contigs via `--contigs` or `--assem_dir`
- Read mapping against assembled contigs with BWA-MEM, filtering for primary alignments with mapping quality ≥ 10
- BAM sorting and indexing with SAMtools
- Optional PCR duplicate marking and removal with Picard MarkDuplicates
- Automatic cleanup of intermediate files and BWA index files

### Output files
- `SAMPLE_sorted.bam`: Sorted BAM file of reads mapped to the assembly (when `--remove_duplicates=f`)
- `SAMPLE_sorted.bam.bai`: Index for the sorted BAM file
- `SAMPLE_sorted_markdup.bam`: Sorted, duplicate-removed BAM file (when `--remove_duplicates=t`)
- `SAMPLE_sorted_markdup.bam.bai`: Index for the duplicate-removed BAM file
- `SAMPLE_sorted_markdup.metrics.txt`: Picard duplicate metrics report (when `--remove_duplicates=t`)
- `SAMPLE.contigs.fa`: Assembled contigs in FASTA format (when MEGAHIT is run)

### Help
```
Usage: 3-assembly_and_map_pipeline.sh [OPTIONS]

Required:
  --reads1 CHAR              Input R1 (or single-end) metagenome reads (fastq/fa)

Paired-end only (ignored when --single_end=t):
  --reads2 CHAR              Input R2 metagenome reads (fastq/fa)

Optional:
  --single_end t|f           Process as single-end reads (default: f)
  --sample_name CHAR         Sample name used to name output files (default: metagenomex)
  --contigs CHAR             Path to pre-assembled contigs file (FASTA format)
                             Supports both compressed (.gz) and uncompressed files
                             Takes precedence over --assem_dir
  --assem_dir CHAR           Directory with previously computed assemblies
                             Will search for: SAMPLE_NAME.contigs.{fa,fasta,fna}[.gz]
                             in ASSEM_DIR/ or ASSEM_DIR/SAMPLE_NAME/
                             Supports both compressed (.gz) and uncompressed files
  --assem_preset CHAR        MEGAHIT preset to generate assembly
                             (default: meta-sensitive)
  --nslots NUM               Number of threads used (default: 12)
  --min_contig_length NUM    Minimum length of contigs to keep (default: 250)
  --output_dir CHAR          Output directory (default: mg-clust_output-1)
  --overwrite t|f            Overwrite previous folder if present (default: f)
  --remove_duplicates t|f    Remove PCR duplicates with Picard (default: f)
  --help                     Print this help and exit
```

### Example usage
```bash
# Paired-end de novo assembly and mapping
./modules/3-assembly_and_map_pipeline.sh \
    --reads1 sample_R1.fastq.gz \
    --reads2 sample_R2.fastq.gz \
    --sample_name sample1 \
    --nslots 16 \
    --output_dir results/sample1_assembly

# Single-end de novo assembly and mapping
./modules/3-assembly_and_map_pipeline.sh \
    --reads1 sample_SE.fastq.gz \
    --single_end t \
    --sample_name sample1 \
    --nslots 16 \
    --output_dir results/sample1_assembly
```

# **Dependencies**

All dependencies are specified in the [environment.yml](environment.yml) file and can be installed via mamba/conda.

**Core tools:**
- [bzip2](http://www.bzip.org) - File compression
- [gzip](https://www.gzip.org) - File compression
- [pigz](https://zlib.net/pigz/) - Parallel gzip compression
- [seqtk](https://github.com/lh3/seqtk) - Sequence processing toolkit
- [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) - Suite including BBDuk and BBMerge for adapter trimming and read merging
- [PEAR](https://cme.h-its.org/exelixis/web/software/pear) - Paired-end read merger
- [fastp](https://github.com/OpenGene/fastp) - Fast all-in-one preprocessing tool
- [MEGAHIT](https://github.com/voutcn/megahit) - De novo assembler for metagenomes
- [BWA](https://github.com/lh3/bwa) - Burrows-Wheeler Aligner for read mapping
- [SAMtools](http://www.htslib.org/) - SAM/BAM file manipulation
- [Picard](https://broadinstitute.github.io/picard/) - Java tools for manipulating sequencing data
- [EMBOSS](http://emboss.sourceforge.net/) - Sequence analysis tools (provides infoseq)

**R and R packages:**
- [R](https://www.r-project.org) - Statistical computing environment
- [tidyverse](https://www.tidyverse.org) - Data manipulation and visualization
- [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) - FASTQ file handling (Bioconductor)
- [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) - Parallel processing
- [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html) - Sequence quality profiling (Bioconductor)
- [optparse](https://cran.r-project.org/web/packages/optparse/index.html) - Command-line argument parsing  

# **License**

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  



