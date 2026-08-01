// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.1: quality-check report (fastp)
// Input:  per-sample reads (R1[, R2])
// Output: per-sample fastp HTML/JSON report + stats + log (diagnostic)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_1_QUALITY_CHECK {

    container "ghcr.io/pereiramemo/mg-proc/1.1-quality-check:${params.container_tag}"
    publishDir "${params.output_dir}/1.1-quality-check-out",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    val single_end_flag
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}"

    script:
    // Nextflow unwraps a single-element path list into a bare scalar (not a List),
    // while 2+ elements become a BlankSeparatedList; normalize back to a List here.
    def rlist = reads instanceof List ? reads : [reads]
    def reads2 = rlist.size() > 1 ? "--reads2 ${rlist[1]}" : ""
    
    """
    1.1-quality-check.py \
        --reads1                    ${rlist[0]} \
        ${reads2} \
        --sample_name               ${sample_name} \
        --single_end                ${single_end_flag ? 't' : 'f'} \
        --output_dir                ${sample_name} \
        --nslots                    ${task.cpus} \
        --min_length                ${params.qc_min_length} \
        --qualified_quality_phred   ${params.qualified_quality_phred} \
        --unqualified_percent_limit ${params.unqualified_percent_limit} \
        --disable_adapter_trimming  ${params.disable_adapter_trimming} \
        --overwrite                 t
    """
}
