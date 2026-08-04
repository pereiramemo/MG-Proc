// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.2: comparative quality-check plots (R / ShortRead / DADA2)
// Input:  all samples' reads staged together (aggregate), plus the input_tsv
//         samplesheet used to map staged files back to sample names
// Output: comparative QC plots + stats + log (diagnostic)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_2_QUALITY_CHECK {

    container "ghcr.io/pereiramemo/mg-proc/1.2-quality-check:${params.container_tag}"
    publishDir "${params.output_dir}/",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean()

    tag "all samples"

    input:
    path reads
    path samplesheet

    output:
    path "1.2-quality-check-out"

    script:
    """
    1.2-quality-check.R \
        --input_dir   . \
        --input_tsv   ${samplesheet} \
        --output_dir  1.2-quality-check-out \
        --sample_size ${params.qc_sample_size} \
        --nslots      ${task.cpus} \
        --overwrite   t
    """
}
