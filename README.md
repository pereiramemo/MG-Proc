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

[1.1-quality_check_fastp.sh](modules/1.1-quality_check_fastp.sh): This script performs comprehensive quality assessment of raw Illumina reads (paired-end or single-end) using fastp. It generates HTML and JSON reports for quality control.

**Main features:**
- Quality assessment without modifying input files
- Paired-end and single-end read support
- HTML and JSON report generation
- Configurable quality and length thresholds
- Batch processing of multiple samples

To see the help run ```./modules/1.1-quality_check_fastp.sh --help```

**Key options:**
- `--input_dir=CHAR`: Directory containing input FASTQ files (required)
- `--output_dir=CHAR`: Directory to output QC reports and plots (required)
- `--single_end=t|f`: Process single-end reads instead of paired-end [default=f]
- `--r1_pattern=CHAR`: Pattern for R1 FASTQ files, or single-end files when `--single_end=t` [default=_R1_001.fastq.gz]
- `--r2_pattern=CHAR`: Pattern for R2 FASTQ files (ignored when `--single_end=t`) [default=_R2_001.fastq.gz]
- `--nslots=NUM`: Number of threads to use [default=12]
- `--min_length=NUM`: Minimum read length filter [default=50]
- `--qualified_quality_phred=NUM`: Minimum quality value [default=20]
- `--overwrite=t|f`: Overwrite previous directory [default=f]

## 1.2-quality_check.R

[1.2-quality_check.R](modules/1.2-quality_check.R): This R script performs comprehensive quality assessment of raw Illumina paired-end reads. It generates multiple plots to evaluate read quality and identify potential issues.

**Main analyses:**
- Calculation of mean quality scores for R1 and R2 reads
- Plotting quality scores versus read counts
- Generating histograms of read count distributions
- Detection and quantification of PhiX contamination

**Output files:**
- `r1_mean_q_vs_nseq.png`: Quality score vs read count for R1
- `r2_mean_q_vs_nseq.png`: Quality score vs read count for R2
- `samples_hist.png`: Histogram of read counts per sample
- `samples_hist_log.png`: Log-transformed histogram of read counts
- `samples_perc_phix_barplot.png`: PhiX contamination levels per sample

To see the help run ```Rscript modules/1.2-quality_check.R --help```

**Key options:**
- `--input_dir=CHARACTER`: Input directory with FASTQ files (required)
- `--output_dir=CHARACTER`: Output directory for plots (required)
- `--nslots=INTEGER`: Number of threads to use [default=12]
- `--r1_pattern=CHARACTER`: Pattern for R1 FASTQ files [default=R1_001.fastq.gz]
- `--r2_pattern=CHARACTER`: Pattern for R2 FASTQ files [default=R2_001.fastq.gz]
- `--overwrite=LOGICAL`: Overwrite previous output [default=FALSE]

## 2-preprocess_pipeline.sh

[2-preprocess_pipeline.sh](modules/2-preprocess_pipeline.sh): This bash pipeline preprocesses raw Illumina reads (paired-end or single-end) from metagenomic samples.

**Main tasks:**
- Optional FASTQ repair to fix malformed records and SRA-format extended headers (using BBTools `reformat.sh` / `repair.sh`)
- Checking for the presence of adapters and trimming them
- Merging paired-end reads using PEAR or BBMerge (PE only)
- Quality trimming of merged, unmerged, and single-end reads
- Computing and plotting read statistics (number and mean length) for intermediate files

**Output files:**
- `*workable.fasta`: FASTA file ready for downstream analyses
- `stats.tsv`: Table with read statistics
- `stats_plots.png`: Plot showing number of sequences and mean read length

Optionally, the pipeline can output quality-checked paired-end reads as FASTQ files.

To see the help run ```./modules/2-preprocess_pipeline.sh --help```

**Required options:**
- `--reads FILE`: R1 (or single-end) FASTQ file
- `--output_dir DIR`: Output directory

**Paired-end only options:**
- `--reads2 FILE`: R2 FASTQ file
- `--merger STR`: pear|bbmerge [default=pear]
- `--min_overlap NUM`: Minimum PE overlap for PEAR [default=10]
- `--output_pe t|f`: Output QC'ed paired-end reads [default=f]
- `--output_merged t|f`: Output merged QC'ed reads [default=t]
- `--pvalue NUM`: p-value for PEAR [default=0.01]

**Optional parameters:**
- `--single_end t|f`: Process as single-end reads [default=f]
- `--repair t|f`: Repair FASTQ before processing [default=f]. For SE: runs `reformat.sh tossbrokenreads=t` to fix malformed records. For PE: first reformats R1/R2 together with `reformat.sh`, then re-pairs with `repair.sh`. Useful for SRA-downloaded data with extended `+` header lines
- `--clean t|f`: Remove intermediates [default=f]
- `--compress t|f`: Compress outputs with pigz [default=f]
- `--min_length NUM`: Minimum read length after trimming [default=75]
- `--min_qual NUM`: Quality trim threshold [default=20]
- `--nslots NUM`: Threads [default=12]
- `--overwrite t|f`: Replace existing output dir [default=f]
- `--plot t|f`: Produce QC plots [default=f]
- `--sample_name STR`: Name prefix [default=metagenomex]
- `--seed NUM`: Random seed for subsampling [default=123]
- `--subsample t|f`: Subsample to 10k reads [default=f]
- `--trim_adapters t|f`: Remove adapters [default=f]

## 3-assembly_and_map_pipeline.sh

[3-assembly_and_map_pipeline.sh](modules/3-assembly_and_map_pipeline.sh): This pipeline performs de novo assembly of metagenomic reads (paired-end or single-end, or using pre-assembled contigs) and maps reads back to the assembled contigs.

**Main tasks:**
- De novo assembly using MEGAHIT (optional — can use pre-assembled contigs). Supports both paired-end (`-1`/`-2`) and single-end (`-r`) MEGAHIT modes
- Read mapping using BWA-MEM with quality filtering (q≥10, primary alignments only)
- BAM file sorting and indexing
- Duplicate removal using Picard (optional)
- Automatic cleanup of intermediate files and BWA indices

**Key features:**
- Paired-end and single-end read support throughout (assembly and mapping)
- Flexible input: run de novo assembly or use existing contigs (supports .gz compressed files)
- Smart contig file discovery with multiple naming patterns
- Quality control: skips processing if < 5 contigs (graceful exit)
- Memory-efficient: direct BAM conversion (avoids large SAM files)
- Clean output: removes intermediate files and BWA indices automatically

To see the help run ```./modules/3-assembly_and_map_pipeline.sh --help```

**Required options:**
- `--reads1 CHAR`: Input R1 (or single-end) metagenome reads (fastq/fa)
- `--sample_name CHAR`: Sample name used to name the files

**Paired-end only options:**
- `--reads2 CHAR`: Input R2 metagenome reads (fastq/fa)

**Optional parameters:**
- `--single_end t|f`: Process as single-end reads [default=f]
- `--contigs CHAR`: Path to pre-assembled contigs file (FASTA format). Supports both compressed (.gz) and uncompressed files. Takes precedence over `--assem_dir`
- `--assem_dir CHAR`: Directory with previously computed assemblies. Will search for `SAMPLE_NAME.contigs.{fa,fasta,fna}[.gz]` in `ASSEM_DIR/` or `ASSEM_DIR/SAMPLE_NAME/`. Supports both compressed (.gz) and uncompressed files
- `--assem_preset CHAR`: MEGAHIT preset to generate assembly [default=meta-sensitive]
- `--nslots NUM`: Number of threads used [default=12]
- `--min_contig_length NUM`: Minimum length of contigs to keep [default=250]
- `--output_dir CHAR`: Output directory [default=mg-clust_output-1]
- `--overwrite t|f`: Overwrite previous folder if present [default=f]
- `--remove_duplicates t|f`: Remove PCR duplicates with Picard [default=f]

**Output files:**
- `SAMPLE_NAME_sorted.bam` (or `SAMPLE_NAME_sorted_markdup.bam` if duplicates removed): Final BAM file
- `SAMPLE_NAME_sorted.bam.bai` (or `SAMPLE_NAME_sorted_markdup.bam.bai`): BAM index
- `SAMPLE_NAME.contigs.fa`: Assembly file (if MEGAHIT was run)
- `SAMPLE_NAME_sorted_markdup.metrics.txt`: Picard duplicate metrics (if duplicates removed)

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



