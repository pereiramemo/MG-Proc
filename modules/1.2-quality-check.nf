// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.2: comparative quality-check plots (R / ShortRead / DADA2)
// Input:  all samples' reads staged together (aggregate)
// Output: comparative QC plots + stats + log (diagnostic)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_2_QUALITY_CHECK {

    container "ghcr.io/pereiramemo/mg-proc/1.2-quality-check:${params.container_tag}"
    publishDir "${params.output_dir}/",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "all samples"

    input:
    path reads

    output:
    path "1.2-quality-check-out"

    script:
    """
    1.2-quality-check.R \
        --input_dir        . \
        --output_dir       1.2-quality-check-out \
        --single_end       ${params.single_end ? 't' : 'f'} \
        --reads_pattern    '${params.reads_pattern}' \
        --se_reads_pattern '${params.se_reads_pattern}' \
        --sample_size      ${params.qc_sample_size} \
        --nslots           ${task.cpus} \
        --overwrite        t
    """
}
