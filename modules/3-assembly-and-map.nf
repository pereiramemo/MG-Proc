// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: de novo assembly + read mapping
// Input:  per-sample QC-trimmed reads from MODULE_2_PREPROCESS
// Output: assembly + sorted (optionally dedup) BAM + stats + log
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_3_ASSEMBLY_AND_MAP {

    container "ghcr.io/pereiramemo/mg-proc/3-assembly-and-map:${params.container_tag}"
    publishDir "${params.output_dir}/3-assembly-and-map-out",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}"

    script:
    def rlist = reads instanceof List ? reads : [reads]
    def reads2    = rlist.size() > 1 ? "--reads2 ${rlist[1]}" : ""
    def contigs   = params.contigs   ? "--contigs ${params.contigs}"     : ""
    def assem_dir = params.assem_dir ? "--assem_dir ${params.assem_dir}" : ""
    """
    3-assembly-and-map.py \
        --reads1            ${rlist[0]} \
        ${reads2} \
        --single_end        ${params.single_end ? 't' : 'f'} \
        --sample_name       ${sample_name} \
        ${contigs} \
        ${assem_dir} \
        --assem_preset      ${params.assem_preset} \
        --min_contig_length ${params.min_contig_length} \
        --remove_duplicates ${params.remove_duplicates} \
        --nslots            ${task.cpus} \
        --output_dir        ${sample_name} \
        --overwrite         t
    """
}
