#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE_1_1_QUALITY_CHECK }   from './modules/1.1-quality-check.nf'
include { MODULE_1_2_QUALITY_CHECK }   from './modules/1.2-quality-check.nf'
include { MODULE_2_PREPROCESS }        from './modules/2-preprocess.nf'
include { MODULE_3_ASSEMBLY_AND_MAP }  from './modules/3-assembly-and-map.nf'

workflow {

    main:
    if (params.help) {
        log.info """
        MG-Proc: metagenomic read processing from raw reads to assembly + mapping

        Usage: nextflow run mg-proc.nf [options]

        General:
          --input_dir         DIR   Input directory with FASTQ files (default: ${params.input_dir})
          --reads_pattern     STR   Paired-end glob for fromFilePairs (default: ${params.reads_pattern})
          --se_reads_pattern  STR   Single-end glob when --single_end true (default: ${params.se_reads_pattern})
          --single_end        BOOL  Process single-end reads (default: ${params.single_end})
          --output_dir        DIR   Output directory (default: ${params.output_dir})
          --nslots            INT   CPU threads per tool (default: ${params.nslots})
          --maxForks          INT   Max parallel process instances (default: ${params.maxForks})
          --full_output       BOOL  Publish all module outputs (default: ${params.full_output})
          --skip_assembly     BOOL  Skip MODULE_3_ASSEMBLY_AND_MAP (default: ${params.skip_assembly})
          --container_tag     STR   Tag of the ghcr.io/pereiramemo/mg-proc/* images (default: ${params.container_tag})

        MODULE_1_1_QUALITY_CHECK — fastp QC report (always runs):
          --qc_min_length             INT  Minimum read length, reporting only (default: ${params.qc_min_length})
          --qualified_quality_phred   INT  Qualified base quality, reporting only (default: ${params.qualified_quality_phred})
          --unqualified_percent_limit INT  Max unqualified base percent (default: ${params.unqualified_percent_limit})
          --disable_adapter_trimming  STR  Disable adapter trimming in report, t/f (default: ${params.disable_adapter_trimming})

        MODULE_1_2_QUALITY_CHECK — comparative QC plots (always runs):
          (finds files via --reads_pattern / --se_reads_pattern from General)
          --qc_sample_size  INT   Reads subsampled per file for QC (default: ${params.qc_sample_size})

        MODULE_2_PREPROCESS — preprocessing:
          --reformat        STR  Reformat FASTQ with reformat.sh, t/f (default: ${params.reformat})
          --repair          STR  Reformat + repair FASTQ (PE), t/f (default: ${params.repair})
          --subsample       STR  Subsample to 10k reads, t/f (default: ${params.subsample})
          --trim_adapters   STR  Remove adapters with BBDuk, t/f (default: ${params.trim_adapters})
          --output_merged   STR  Merge PE reads, t/f (default: ${params.output_merged})
          --merger          STR  pear | bbmerge (default: ${params.merger})
          --min_overlap     INT  Minimum PE overlap for PEAR (default: ${params.min_overlap})
          --pvalue          NUM  p-value for PEAR (default: ${params.pvalue})
          --min_length      INT  Minimum read length after trimming (default: ${params.min_length})
          --min_qual        INT  Quality trim threshold (default: ${params.min_qual})
          --seed            INT  Random seed for subsampling (default: ${params.seed})
          --clean           STR  Remove intermediates, t/f (default: ${params.clean})
          --compress        STR  Compress outputs with pigz, t/f (default: ${params.compress})

        MODULE_3_ASSEMBLY_AND_MAP — assembly + mapping:
          --assem_preset      STR  MEGAHIT preset (default: ${params.assem_preset})
          --min_contig_length INT  Minimum contig length to keep (default: ${params.min_contig_length})
          --remove_duplicates STR  Remove PCR duplicates with Picard, t/f (default: ${params.remove_duplicates})
          --contigs           STR  Pre-assembled contigs FASTA (default: none)
          --assem_dir         STR  Directory with previous assemblies (default: none)
        """.stripIndent()
        exit 0
    }

    // Build the reads channel: paired-end (fromFilePairs) or single-end (fromPath).
    if (params.single_end) {
        reads_ch = channel.fromPath(
            "${params.input_dir}/${params.se_reads_pattern}", checkIfExists: true
        ).map { f -> tuple(f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), [f]) }
    } else {
        reads_ch = channel.fromFilePairs(
            "${params.input_dir}/${params.reads_pattern}", checkIfExists: true
        )
    }

    log.info "Read type: ${params.single_end ? 'single-end' : 'paired-end'}"
    if (params.skip_assembly) {
        log.info "Assembly + mapping will be skipped"
    }

    // MODULE_1_1_QUALITY_CHECK: per-sample fastp report (diagnostic, always runs)
    MODULE_1_1_QUALITY_CHECK(reads_ch)

    // MODULE_1_2_QUALITY_CHECK: comparative plots over all samples (diagnostic, always runs)
    all_reads = reads_ch.map { _sample_name, r -> r }.flatten().collect()
    MODULE_1_2_QUALITY_CHECK(all_reads)

    // MODULE_2_PREPROCESS: per-sample preprocessing
    preprocess = MODULE_2_PREPROCESS(reads_ch)

    // MODULE_3_ASSEMBLY_AND_MAP: assemble + map the preprocessed QC-trimmed reads
    if (!params.skip_assembly) {
        qc_reads = preprocess.qc_reads.map { sample_name, r ->
            tuple(sample_name, r instanceof List ? r : [r])
        }
        MODULE_3_ASSEMBLY_AND_MAP(qc_reads)
    }
}
