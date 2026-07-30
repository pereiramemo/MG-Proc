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
          --input_tsv         FILE  TSV samplesheet: sample_name, reads1, reads2
                                     (empty reads2 => single-end; default: ${params.input_tsv})
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
          --qc_sample_size  INT   Reads subsampled per file for QC (default: ${params.qc_sample_size})

        MODULE_2_PREPROCESS — preprocessing:
          --reformat        STR  Reformat FASTQ with reformat.sh, t/f (default: ${params.reformat})
          --repair          STR  Reformat + repair FASTQ (PE), t/f (default: ${params.repair})
          --subsample       STR  Subsample to 10k reads, t/f (default: ${params.subsample})
          --trim_adapters   STR  Remove adapters with BBDuk, t/f (default: ${params.trim_adapters})
          --output_pe       STR  Output PE R1/R2 QC reads; only takes effect with
                                  --skip_assembly true, t/f (default: ${params.output_pe})
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

    // Build the reads channel from the TSV samplesheet (sample_name, reads1,
    // reads2). An empty reads2 marks a single-end sample; the sheet must be
    // all-PE or all-SE (mixed sheets are rejected below). Parsed eagerly
    // (file().splitCsv(), not a channel) since single_end has to be known
    // before any process is invoked.
    samplesheet = file(params.input_tsv)
    rows = samplesheet.splitCsv(header: true, sep: '\t')
    if (rows.isEmpty()) {
        error "input_tsv '${params.input_tsv}' has no data rows"
    }

    single_end_flags_list = rows.collect { (it.reads2?.trim()) ? false : true }.unique()
    if (single_end_flags_list.size() > 1) {
        error "input_tsv mixes single-end and paired-end rows (reads2 must be either always empty or always populated)"
    }
    single_end_flag = single_end_flags_list[0]

    reads_ch = channel.fromList(rows).map { row ->
        // def scopes files to this closure call; without it Groovy would bind
        // files in the enclosing script scope, shared/reused across every row.
        def files = single_end_flag ? [file(row.reads1)] : [file(row.reads1), file(row.reads2)]
        tuple(row.sample_name, files)
    }

    log.info "Read type: ${single_end_flag ? 'single-end' : 'paired-end'}"
    if (params.skip_assembly) {
        log.info "Assembly + mapping will be skipped"
    }

    // MODULE_1_1_QUALITY_CHECK: per-sample fastp report (diagnostic, always runs)
    MODULE_1_1_QUALITY_CHECK(single_end_flag, reads_ch)

    // MODULE_1_2_QUALITY_CHECK: comparative plots over all samples (diagnostic, always runs)
    all_reads = reads_ch.map { _sample_name, r -> r }.flatten().collect()
    MODULE_1_2_QUALITY_CHECK(all_reads, samplesheet)

    // MODULE_2_PREPROCESS: per-sample preprocessing
    preprocess = MODULE_2_PREPROCESS(single_end_flag, reads_ch)

    // MODULE_3_ASSEMBLY_AND_MAP: assemble + map the preprocessed QC-trimmed reads
    if (!params.skip_assembly) {
        qc_reads = preprocess.qc_reads.map { sample_name, r ->
            tuple(sample_name, r instanceof List ? r : [r])
        }
        MODULE_3_ASSEMBLY_AND_MAP(single_end_flag, qc_reads)
    }
}
