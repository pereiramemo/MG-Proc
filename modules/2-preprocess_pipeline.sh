#!/usr/bin/env bash
###############################################################################
# preprocess_pipeline.sh
###############################################################################

set -euo pipefail

###############################################################################
# 1. Environment
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/conf.sh"

###############################################################################
# 2. Help
###############################################################################

show_usage() {
cat <<EOF
Usage: $(basename "$0") [OPTIONS]

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
  --repair t|f          Repair FASTQ before processing (default f)
                        SE: reformat.sh fixes malformed records
                        PE: repair.sh re-pairs mismatched reads
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
EOF
}

###############################################################################
# 3. Parse arguments with getopt
###############################################################################

CLEAN="f"
COMPRESS="f"
MERGER="pear"
MIN_LENGTH=75
MIN_OVERLAP=10
MIN_QUAL=20
NSLOTS=12
OUTPUT_DIR=""
OUTPUT_PE="f"
OUTPUT_MERGED="t"
OVERWRITE="f"
PVALUE="0.01"
PLOT="f"
R1=""
R2=""
SAMPLE_NAME="metagenomex"
SEED=123
SINGLE_END="f"
REPAIR="f"
SUBSAMPLE="f"
TRIM_ADAPTERS="f"

# Initialize optional path vars to avoid -u issues later
R1_UNCOMPRESSED=""
R2_UNCOMPRESSED=""
R1_REDU=""
R2_REDU=""
R1_AT=""
R2_AT=""
R1_QC=""
R2_QC=""
R_ASSEM=""
R1_UNASSEM=""
R2_UNASSEM=""
R_DISCARD=""
R_ASSEM_QC=""
R1_UNASSEM_QC=""
R2_UNASSEM_QC=""
R_ASSEM_QC_FA=""
R1_UNASSEM_QC_FA=""
R2_UNASSEM_QC_FA=""
SE_QC=""
SE_QC_FA=""
R1_REFORMATTED=""
R2_REFORMATTED=""
R1_REPAIRED=""
R2_REPAIRED=""
SINGLETON_FILE=""

check_cmd getopt

ARGS=$(getopt -o '' \
  --long help,clean:,compress:,merger:,min_length:,min_overlap:,min_qual:,nslots:,\
output_dir:,output_pe:,output_merged:,overwrite:,pvalue:,plot:,reads:,reads2:,\
sample_name:,seed:,single_end:,repair:,subsample:,trim_adapters: \
  -n "$(basename "$0")" -- "$@" \
  ) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    --help) show_usage; exit 0 ;;
    --clean) CLEAN="$2"; shift 2 ;;
    --compress) COMPRESS="$2"; shift 2 ;;
    --merger) MERGER="$2"; shift 2 ;;
    --min_length) MIN_LENGTH="$2"; shift 2 ;;
    --min_overlap) MIN_OVERLAP="$2"; shift 2 ;;
    --min_qual) MIN_QUAL="$2"; shift 2 ;;
    --nslots) NSLOTS="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --output_pe) OUTPUT_PE="$2"; shift 2 ;;
    --output_merged) OUTPUT_MERGED="$2"; shift 2 ;;
    --overwrite) OVERWRITE="$2"; shift 2 ;;
    --pvalue) PVALUE="$2"; shift 2 ;;
    --plot) PLOT="$2"; shift 2 ;;
    --reads) R1="$2"; shift 2 ;;
    --reads2) R2="$2"; shift 2 ;;
    --repair) REPAIR="$2"; shift 2 ;;
    --sample_name) SAMPLE_NAME="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --single_end) SINGLE_END="$2"; shift 2 ;;
    --subsample) SUBSAMPLE="$2"; shift 2 ;;
    --trim_adapters) TRIM_ADAPTERS="$2"; shift 2 ;;
    --) shift; break ;;
    *) log_error "Internal getopt error"; exit 1 ;;
  esac
done

###############################################################################
# 4. Validate arguments
###############################################################################

if [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "--output_dir is required."
    exit 1
fi

if [[ -z "${R1}" ]]; then
    log_error "--reads is required."
    exit 1
fi

check_file "${R1}"

if [[ "${SINGLE_END}" == "f" ]]; then
    if [[ -z "${R2}" ]]; then
        log_error "--reads2 is required for paired-end mode."
        exit 1
    fi
    check_file "${R2}"
fi

# Validate boolean flags
for flag in CLEAN COMPRESS OUTPUT_PE OUTPUT_MERGED OVERWRITE SINGLE_END REPAIR SUBSAMPLE TRIM_ADAPTERS PLOT; do
    if ! [[ "${!flag}" =~ ^[tf]$ ]]; then
        log_error "Flag --$(echo "${flag}" | tr 'A-Z' 'a-z') must be 't' or 'f' (got '${!flag}')."
        exit 1
    fi
done

# Validate numeric parameters
for num in MIN_LENGTH MIN_OVERLAP MIN_QUAL NSLOTS SEED; do
    if ! [[ "${!num}" =~ ^[0-9]+$ ]]; then
        log_error "Parameter ${num} must be numeric (got '${!num}')."
        exit 1
    fi
done

# Validate PVALUE (float)
if ! [[ "${PVALUE}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    log_error "Parameter PVALUE must be numeric (got '${PVALUE}')."
    exit 1
fi

# Validate merger (PE only)
if [[ "${SINGLE_END}" == "f" ]]; then
    if [[ "${MERGER}" != "pear" && "${MERGER}" != "bbmerge" ]]; then
        log_error "--merger must be 'pear' or 'bbmerge' (got '${MERGER}')."
        exit 1
    fi
fi

###############################################################################
# 5. Check dependencies
###############################################################################

check_dependencies

###############################################################################
# 6. Prepare output directory
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ "${OVERWRITE}" == "t" ]]; then
        log_warn "Overwriting existing directory: ${OUTPUT_DIR}"
        rm -rf "${OUTPUT_DIR}"
    else
        log_error "Output directory exists: ${OUTPUT_DIR}"
        exit 1
    fi
fi

mkdir -p "${OUTPUT_DIR}"

###############################################################################
# 7. Detect compression and prepare inputs
###############################################################################

log "Checking input compression..."
UNCOMPRESSED="f"

TEST_FILE_BZIP2=$(file "${R1}" | grep -E "bzip2" || true)
TEST_FILE_GZIP=$(file "${R1}" | grep -E "gzip" || true)

if [[ -n "${TEST_FILE_BZIP2}" ]]; then
    log "Uncompressing bzip2 ..."
    R1_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
    if ! "${bzip2}" --decompress --keep --stdout "${R1}" > "${R1_UNCOMPRESSED}"; then
        log_error "bzip2 R1 failed"
        exit 1
    fi
    R1="${R1_UNCOMPRESSED}"

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
        if ! "${bzip2}" --decompress --keep --stdout "${R2}" > "${R2_UNCOMPRESSED}"; then
            log_error "bzip2 R2 failed"
            exit 1
        fi
        R2="${R2_UNCOMPRESSED}"
    fi

    UNCOMPRESSED="t"
elif [[ -n "${TEST_FILE_GZIP}" ]]; then
    log "Uncompressing gzip ..."
    R1_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
    if ! "${gunzip}" --keep --stdout "${R1}" > "${R1_UNCOMPRESSED}"; then
        log_error "gunzip R1 failed"
        exit 1
    fi
    R1="${R1_UNCOMPRESSED}"

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
        if ! "${gunzip}" --keep --stdout "${R2}" > "${R2_UNCOMPRESSED}"; then
            log_error "gunzip R2 failed"
            exit 1
        fi
        R2="${R2_UNCOMPRESSED}"
    fi

    UNCOMPRESSED="t"
else
    log "Linking original FASTQ files..."
    R1_ABS="$(realpath "${R1}")"

    if ! ln -s "${R1_ABS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"; then
        log_error "ln -s R1 failed"
        exit 1
    fi
    R1="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2_ABS="$(realpath "${R2}")"
        if ! ln -s "${R2_ABS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"; then
            log_error "ln -s R2 failed"
            exit 1
        fi
        R2="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
    fi
fi

###############################################################################
# 8. Repair FASTQ files (optional)
###############################################################################

if [[ "${REPAIR}" == "t" ]]; then
    if [[ "${SINGLE_END}" == "t" ]]; then
        log "Reformatting single-end FASTQ with reformat.sh..."
        R1_REFORMATTED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_reformatted-00.fastq"

        if ! "${reformat}" \
            in="${R1}" \
            out="${R1_REFORMATTED}" \
            tossbrokenreads=t; then
            log_error "reformat.sh reformat of R1 failed"
            exit 1
        fi
        R1="${R1_REFORMATTED}"
    else
        log "Reformatting paired-end FASTQ with reformat.sh..."
        R1_REFORMATTED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_reformatted-00.fastq"
        R2_REFORMATTED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_reformatted-00.fastq"

        if ! "${reformat}" \
            in="${R1}" in2="${R2}" \
            out="${R1_REFORMATTED}" out2="${R2_REFORMATTED}" \
            tossbrokenreads=t; then
            log_error "reformat.sh reformat of R1/R2 failed"
            exit 1
        fi

        log "Repairing paired-end FASTQ with repair.sh..."
        R1_REPAIRED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_repaired-00.fastq"
        R2_REPAIRED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_repaired-00.fastq"
        SINGLETON_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_singletons_repaired-00.fastq"

        if ! "${repair_sh}" \
            in="${R1_REFORMATTED}" in2="${R2_REFORMATTED}" \
            out="${R1_REPAIRED}" out2="${R2_REPAIRED}" \
            outs="${SINGLETON_FILE}"; then
            log_error "repair.sh failed"
            exit 1
        fi
        R1="${R1_REPAIRED}"
        R2="${R2_REPAIRED}"
    fi
fi

###############################################################################
# 9. Subsample data (optional)
###############################################################################

if [[ "${SUBSAMPLE}" == "t" ]]; then
    log "Subsampling to 10,000 reads (seed: ${SEED}) ..."
    R1_REDU="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_redu-00.fastq"

    if ! "${seqtk}" sample -s"${SEED}" "${R1}" 10000 > "${R1_REDU}"; then
        log_error "seqtk subsampling R1 failed"
        exit 1
    fi
    R1="${R1_REDU}"

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2_REDU="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_redu-00.fastq"
        if ! "${seqtk}" sample -s"${SEED}" "${R2}" 10000 > "${R2_REDU}"; then
            log_error "seqtk subsampling R2 failed"
            exit 1
        fi
        R2="${R2_REDU}"
    fi
fi

###############################################################################
# 10. Adapter trimming (optional)
###############################################################################

if [[ "${TRIM_ADAPTERS}" == "t" ]]; then
    log "Trimming adapters ..."
    R1_AT="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_at-01.fastq"

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2_AT="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_at-01.fastq"
        if ! "${bbduk}" \
            in="${R1}" in2="${R2}" \
            out="${R1_AT}" out2="${R2_AT}" \
            threads="${NSLOTS}" \
            ktrim=r k=23 mink=11 hdist=1 tpe tbo \
            ref="${ADAPTERS}"; then
            log_error "bbduk adapter trimming failed"
            exit 1
        fi
        R2="${R2_AT}"
    else
        if ! "${bbduk}" \
            in="${R1}" \
            out="${R1_AT}" \
            threads="${NSLOTS}" \
            ktrim=r k=23 mink=11 hdist=1 \
            ref="${ADAPTERS}"; then
            log_error "bbduk adapter trimming failed"
            exit 1
        fi
    fi

    R1="${R1_AT}"
fi

###############################################################################
# 11. Quality trim paired-end reads (optional, PE only)
###############################################################################

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_PE}" == "t" ]]; then
    log "Quality trimming paired-end reads ..."
    R1_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_qc-02.fastq"
    R2_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_qc-02.fastq"

    if ! "${bbduk}" \
        in="${R1}" in2="${R2}" \
        out="${R1_QC}" out2="${R2_QC}" \
        minlength="${MIN_LENGTH}" \
        threads="${NSLOTS}" \
        qtrim=rl trimq="${MIN_QUAL}"; then
        log_error "bbduk quality trimming paired-end reads failed"
        exit 1
    fi

    if [[ "${COMPRESS}" == "t" ]]; then
        if ! "${pigz}" --keep --processes "${NSLOTS}" "${R1_QC}" "${R2_QC}"; then
            log_error "pigz compressing ${R1_QC} and ${R2_QC} failed"
            exit 1
        fi
    fi
fi

###############################################################################
# 12. Merge paired-end reads (optional, PE only)
###############################################################################

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" ]]; then
    if [[ "${MERGER}" == "pear" ]]; then
        log "Merging reads with PEAR..."
        if ! "${pear}" \
            -f "${R1}" -r "${R2}" \
            -o "${OUTPUT_DIR}/${SAMPLE_NAME}" \
            -j "${NSLOTS}" \
            -p "${PVALUE}" \
            --min-overlap "${MIN_OVERLAP}"; then
            log_error "Merge with ${pear} failed"
            exit 1
        fi

        mv "${OUTPUT_DIR}/${SAMPLE_NAME}.assembled.fastq" \
           "${OUTPUT_DIR}/${SAMPLE_NAME}_assembled-02.fastq"
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}.unassembled.forward.fastq" \
           "${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R1-02.fastq"
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}.unassembled.reverse.fastq" \
           "${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R2-02.fastq"
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}.discarded.fastq" \
           "${OUTPUT_DIR}/${SAMPLE_NAME}_discarded-02.fastq"

        R_ASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_assembled-02.fastq"
        R1_UNASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R1-02.fastq"
        R2_UNASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R2-02.fastq"
        R_DISCARD="${OUTPUT_DIR}/${SAMPLE_NAME}_discarded-02.fastq"

        if [[ ! -s "${R_ASSEM}" ]]; then
            log_error "Failed merging: 0 merged reads with ${pear}"
            exit 1
        fi
    elif [[ "${MERGER}" == "bbmerge" ]]; then
        log "Merging reads with BBMerge..."
        R_ASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_assembled-02.fastq"
        R1_UNASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R1-02.fastq"
        R2_UNASSEM="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R2-02.fastq"
        INSERT_HIST="${OUTPUT_DIR}/${SAMPLE_NAME}_insert_hist.txt"
        R_DISCARD="${OUTPUT_DIR}/${SAMPLE_NAME}_discarded-02.fastq"

        if ! "${bbmerge}" \
            in="${R1}" in2="${R2}" \
            out="${R_ASSEM}" \
            outu="${R1_UNASSEM}" outu2="${R2_UNASSEM}" \
            ihist="${INSERT_HIST}"; then
            log_error "Merge with ${bbmerge} failed"
            exit 1
        fi

        touch "${R_DISCARD}"
    fi
fi

###############################################################################
# 13. Quality check merged reads (PE only)
###############################################################################

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" && -s "${R_ASSEM}" ]]; then
    log "Quality trimming merged reads..."
    R_ASSEM_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_assembled_qc-03.fastq"

    if ! "${bbduk}" \
        in="${R_ASSEM}" \
        out="${R_ASSEM_QC}" \
        minlength="${MIN_LENGTH}" \
        threads="${NSLOTS}" \
        qtrim=rl trimq="${MIN_QUAL}"; then
        log_error "bbduk quality trimming assembled reads failed"
        exit 1
    fi
fi

###############################################################################
# 14. Quality check unmerged reads (PE only)
###############################################################################

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" && -s "${R1_UNASSEM}" ]]; then
    log "Quality trimming unmerged reads..."
    R1_UNASSEM_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R1_qc-03.fastq"
    R2_UNASSEM_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R2_qc-03.fastq"

    if ! "${bbduk}" \
        in="${R1_UNASSEM}" in2="${R2_UNASSEM}" \
        out="${R1_UNASSEM_QC}" out2="${R2_UNASSEM_QC}" \
        minlength="${MIN_LENGTH}" \
        threads="${NSLOTS}" \
        qtrim=rl trimq="${MIN_QUAL}"; then
        log_error "bbduk quality trimming unassembled reads failed"
        exit 1
    fi
fi

###############################################################################
# 15. Quality trim single-end reads (SE only)
###############################################################################

if [[ "${SINGLE_END}" == "t" ]]; then
    log "Quality trimming single-end reads..."
    SE_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_se_qc-02.fastq"

    if ! "${bbduk}" \
        in="${R1}" \
        out="${SE_QC}" \
        minlength="${MIN_LENGTH}" \
        threads="${NSLOTS}" \
        qtrim=rl trimq="${MIN_QUAL}"; then
        log_error "bbduk quality trimming single-end reads failed"
        exit 1
    fi

    if [[ "${COMPRESS}" == "t" ]]; then
        if ! "${pigz}" --keep --processes "${NSLOTS}" "${SE_QC}"; then
            log_error "pigz compressing ${SE_QC} failed"
            exit 1
        fi
    fi
fi

###############################################################################
# 16. Convert to FASTA
###############################################################################

log "Converting to FASTA format..."

if [[ "${SINGLE_END}" == "t" && -s "${SE_QC}" ]]; then
    SE_QC_FA="${OUTPUT_DIR}/${SAMPLE_NAME}_se_qc-02.fasta"

    if ! "${fq2fa}" "${SE_QC}" > "${SE_QC_FA}"; then
        log_error "fq2fa single-end reads failed"
        exit 1
    fi

    if [[ "${COMPRESS}" == "t" ]]; then
        if ! "${pigz}" --keep --processes "${NSLOTS}" "${SE_QC_FA}"; then
            log_error "pigz compressing ${SE_QC_FA} failed"
            exit 1
        fi
    fi
fi

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" && -s "${R_ASSEM_QC}" ]]; then
    R_ASSEM_QC_FA="${OUTPUT_DIR}/${SAMPLE_NAME}_assembled_qc-03.fasta"

    if ! "${fq2fa}" "${R_ASSEM_QC}" > "${R_ASSEM_QC_FA}"; then
        log_error "fq2fa merged reads failed"
        exit 1
    fi

    if [[ "${COMPRESS}" == "t" ]]; then
        if ! "${pigz}" --keep --processes "${NSLOTS}" "${R_ASSEM_QC_FA}"; then
            log_error "pigz compressing ${R_ASSEM_QC_FA} failed"
            exit 1
        fi
    fi
fi

if [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" && -s "${R1_UNASSEM_QC}" ]]; then
    R1_UNASSEM_QC_FA="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R1_qc-03.fasta"
    R2_UNASSEM_QC_FA="${OUTPUT_DIR}/${SAMPLE_NAME}_unassembled_R2_qc-03.fasta"

    if ! "${fq2fa}" "${R1_UNASSEM_QC}" > "${R1_UNASSEM_QC_FA}"; then
        log_error "fq2fa unmerged forward reads failed"
        exit 1
    fi
    if ! "${fq2fa}" "${R2_UNASSEM_QC}" > "${R2_UNASSEM_QC_FA}"; then
        log_error "fq2fa unmerged reverse reads failed"
        exit 1
    fi

    if [[ "${COMPRESS}" == "t" ]]; then
        if ! "${pigz}" --keep --processes "${NSLOTS}" "${R1_UNASSEM_QC_FA}" "${R2_UNASSEM_QC_FA}"; then
            log_error "pigz compressing ${R1_UNASSEM_QC_FA} and ${R2_UNASSEM_QC_FA} failed"
            exit 1
        fi
    fi
fi

###############################################################################
# 17. Compute stats
###############################################################################

log "Computing statistics ..."

STATS_TMP="${OUTPUT_DIR}/stats_tmp.tsv"
> "${STATS_TMP}"

while IFS= read -r FILE; do
    [[ -z "${FILE}" ]] && continue

    N=$(count_fastq "${FILE}")
    if [[ -z "${N}" ]]; then
        log_error "count_fastq failed for ${FILE}"
        exit 1
    fi

    if [[ -s "${FILE}" ]]; then
        L=$(compute_mean_length "${FILE}")
        if [[ -z "${L}" ]]; then
            log_error "compute_mean_length failed for ${FILE}"
            exit 1
        fi
    else
        L=0
    fi

    BASENAME="$(basename "${FILE}")"
    echo -e "${SAMPLE_NAME}\t${BASENAME}\tnum_seq\t${N}" >> "${STATS_TMP}"
    echo -e "${SAMPLE_NAME}\t${BASENAME}\tmean_length\t${L}" >> "${STATS_TMP}"
done < <(find "${OUTPUT_DIR}" -name "*.fastq")

awk -v OFS="\t" '{
    n = gensub(".*-([0-9]+).fastq","\\1","g",$2);
    print $0, n
}' "${STATS_TMP}" > "${OUTPUT_DIR}/stats.tsv"

rm -f "${STATS_TMP}"

###############################################################################
# 18. Plot stats (optional)
###############################################################################

if [[ "${PLOT}" == "t" ]]; then
    log "Creating plot ..."
    if ! "${rscript}" --vanilla \
        "${plots}" \
        "${OUTPUT_DIR}/stats.tsv" \
        "${OUTPUT_DIR}/stats_plots.png"; then
        log_error "R plotting failed"
        exit 1
    fi
fi

###############################################################################
# 19. Clean intermediates (optional)
###############################################################################

if [[ "${CLEAN}" == "t" ]]; then
    log "Cleaning intermediate files ..."

    # 00 files
    if [[ "${UNCOMPRESSED}" == "t" ]]; then
        [[ -n "${R1_UNCOMPRESSED}" && -e "${R1_UNCOMPRESSED}" ]] && rm -f "${R1_UNCOMPRESSED}"
        [[ -n "${R2_UNCOMPRESSED}" && -e "${R2_UNCOMPRESSED}" ]] && rm -f "${R2_UNCOMPRESSED}"
    fi

    if [[ "${SUBSAMPLE}" == "t" ]]; then
        [[ -n "${R1_REDU}" && -e "${R1_REDU}" ]] && rm -f "${R1_REDU}"
        [[ -n "${R2_REDU}" && -e "${R2_REDU}" ]] && rm -f "${R2_REDU}"
    fi

    if [[ "${REPAIR}" == "t" ]]; then
        [[ -n "${R1_REFORMATTED}" && -e "${R1_REFORMATTED}" ]] && rm -f "${R1_REFORMATTED}"
        [[ -n "${R2_REFORMATTED}" && -e "${R2_REFORMATTED}" ]] && rm -f "${R2_REFORMATTED}"
        [[ -n "${R1_REPAIRED}" && -e "${R1_REPAIRED}" ]] && rm -f "${R1_REPAIRED}"
        [[ -n "${R2_REPAIRED}" && -e "${R2_REPAIRED}" ]] && rm -f "${R2_REPAIRED}"
        [[ -n "${SINGLETON_FILE}" && -e "${SINGLETON_FILE}" ]] && rm -f "${SINGLETON_FILE}"
    fi

    # 01 files
    if [[ "${TRIM_ADAPTERS}" == "t" ]]; then
        [[ -n "${R1_AT}" && -e "${R1_AT}" ]] && rm -f "${R1_AT}"
        [[ -n "${R2_AT}" && -e "${R2_AT}" ]] && rm -f "${R2_AT}"
    fi

    # 02 files
    if [[ "${SINGLE_END}" == "t" ]]; then
        if [[ "${COMPRESS}" == "t" ]]; then
            [[ -n "${SE_QC}" && -f "${SE_QC}" ]] && rm -f "${SE_QC}"
        fi
    else
        if [[ "${OUTPUT_PE}" == "t" && "${COMPRESS}" == "t" ]]; then
            [[ -n "${R1_QC}" && -f "${R1_QC}" ]] && rm -f "${R1_QC}"
            [[ -n "${R2_QC}" && -f "${R2_QC}" ]] && rm -f "${R2_QC}"
        fi

        if [[ "${OUTPUT_MERGED}" == "t" ]]; then
            for f in "${R_ASSEM}" "${R1_UNASSEM}" "${R2_UNASSEM}" "${R_DISCARD}" \
                     "${R_ASSEM_QC}" "${R1_UNASSEM_QC}" "${R2_UNASSEM_QC}"; do
                [[ -n "${f}" && -e "${f}" ]] && rm -f "${f}"
            done
            
            # 03 files
            if [[ "${COMPRESS}" == "t" ]]; then
                for f in "${R1_UNASSEM_QC_FA}" "${R2_UNASSEM_QC_FA}" "${R_ASSEM_QC_FA}"; do
                    [[ -n "${f}" && -e "${f}" ]] && rm -f "${f}"
                done
            fi
        fi
    fi
fi

###############################################################################
# 20. Rename to workable
###############################################################################

if [[ "${SINGLE_END}" == "t" && -n "${SE_QC_FA}" ]]; then
    if [[ "${COMPRESS}" == "t" ]]; then
        if ! mv "${SE_QC_FA}.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_workable.fasta.gz"; then
            log_error "Rename to workable (compressed) failed"
            exit 1
        fi
    else
        if ! mv "${SE_QC_FA}" "${OUTPUT_DIR}/${SAMPLE_NAME}_workable.fasta"; then
            log_error "Rename to workable (uncompressed) failed"
            exit 1
        fi
    fi
elif [[ "${SINGLE_END}" == "f" && "${OUTPUT_MERGED}" == "t" && -n "${R_ASSEM_QC_FA}" ]]; then
    if [[ "${COMPRESS}" == "t" ]]; then
        if ! mv "${R_ASSEM_QC_FA}.gz" "${OUTPUT_DIR}/${SAMPLE_NAME}_workable.fasta.gz"; then
            log_error "Rename to workable (compressed) failed"
            exit 1
        fi
    else
        if ! mv "${R_ASSEM_QC_FA}" "${OUTPUT_DIR}/${SAMPLE_NAME}_workable.fasta"; then
            log_error "Rename to workable (uncompressed) failed"
            exit 1
        fi
    fi
fi

###############################################################################
# 21. End
###############################################################################

log "${GREEN}preprocess_pipeline.sh exited successfully${NC}"
