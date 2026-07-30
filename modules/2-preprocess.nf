// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: preprocessing pipeline
// Input:  per-sample reads (R1[, R2])
// Output: workable FASTA + stats/log; emits QC-trimmed reads for assembly
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_PREPROCESS {

    container "ghcr.io/pereiramemo/mg-proc/2-preprocess:${params.container_tag}"
    publishDir "${params.output_dir}/2-preprocess-out",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    val single_end_flag
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}/output/*_qc-02.fastq*"), emit: qc_reads, optional: true
    path "${sample_name}",                                                emit: dir

    script:
    // Nextflow unwraps a single-element path list into a bare scalar (not a List),
    // while 2+ elements become a BlankSeparatedList; normalize back to a List here.
    def rlist = reads instanceof List ? reads : [reads]
    def reads2 = rlist.size() > 1 ? "--reads2 ${rlist[1]}" : ""
    // Force PE QC reads on only when MODULE_3_ASSEMBLY_AND_MAP will consume them;
    // otherwise honor params.output_pe. SE always emits _se_qc-02.fastq.
    def output_pe = single_end_flag ? 'f' : (params.skip_assembly ? params.output_pe : 't')
    """
    2-preprocess.py \
        --reads         ${rlist[0]} \
        ${reads2} \
        --sample_name   ${sample_name} \
        --single_end    ${single_end_flag ? 't' : 'f'} \
        --output_dir    ${sample_name} \
        --reformat      ${params.reformat} \
        --repair        ${params.repair} \
        --subsample     ${params.subsample} \
        --trim_adapters ${params.trim_adapters} \
        --output_pe     ${output_pe} \
        --output_merged ${params.output_merged} \
        --merger        ${params.merger} \
        --min_overlap   ${params.min_overlap} \
        --pvalue        ${params.pvalue} \
        --min_length    ${params.min_length} \
        --min_qual      ${params.min_qual} \
        --seed          ${params.seed} \
        --clean         ${params.clean} \
        --compress      ${params.compress} \
        --nslots        ${task.cpus} \
        --overwrite     t
    """
}
