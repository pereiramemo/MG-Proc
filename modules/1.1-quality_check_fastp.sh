#!/usr/bin/env bash

###############################################################################
# quality_check_fastp.sh
###############################################################################

set -euo pipefail

###############################################################################
# 1. Environment
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/conf.sh"

###############################################################################
# 2. Define help
###############################################################################

show_usage() {
    cat << EOF
Usage: $(basename "$0") <options>

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
EOF
}

###############################################################################
# 3. Define default parameters
###############################################################################

INPUT_DIR=""
OUTPUT_DIR=""
SINGLE_END="f"
R1_PATTERN="_R1_001.fastq.gz"
R2_PATTERN="_R2_001.fastq.gz"
NSLOTS=12
MIN_LENGTH=50
QUALIFIED_QUALITY_PHRED=20
UNQUALIFIED_PERCENT_LIMIT=40
DISABLE_ADAPTER_TRIMMING="t"
HTML_REPORT="t"
JSON_REPORT="t"
OVERWRITE="f"

###############################################################################
# 4. Parse arguments
###############################################################################

if [[ $# -eq 0 ]]; then
    show_usage
    exit 1
fi

ARGS=$(
  getopt -o h \
  -l help,input_dir:,output_dir:,single_end:,r1_pattern:,r2_pattern:,nslots:,min_length:,\
qualified_quality_phred:,unqualified_percent_limit:,disable_adapter_trimming:,\
html_report:,json_report:,overwrite: \
  -n "$(basename "$0")" -- "$@" \
  ) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    -h|--help) show_usage; exit 0 ;;
    --input_dir) INPUT_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --single_end) SINGLE_END="$2"; shift 2 ;;
    --r1_pattern) R1_PATTERN="$2"; shift 2 ;;
    --r2_pattern) R2_PATTERN="$2"; shift 2 ;;
    --nslots) NSLOTS="$2"; shift 2 ;;
    --min_length) MIN_LENGTH="$2"; shift 2 ;;
    --qualified_quality_phred) QUALIFIED_QUALITY_PHRED="$2"; shift 2 ;;
    --unqualified_percent_limit) UNQUALIFIED_PERCENT_LIMIT="$2"; shift 2 ;;
    --disable_adapter_trimming) DISABLE_ADAPTER_TRIMMING="$2"; shift 2 ;;
    --html_report) HTML_REPORT="$2"; shift 2 ;;
    --json_report) JSON_REPORT="$2"; shift 2 ;;
    --overwrite) OVERWRITE="$2"; shift 2 ;;
    --) shift; break ;;
    *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
  esac
done

###############################################################################
# 5. Validate parameters
###############################################################################

# Required parameters
if [[ -z "${INPUT_DIR}" ]]; then
    log_error "--input_dir is required."
    exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
    log_error "--output_dir is required."
    exit 1
fi

# Check input directory exists
if [[ ! -d "${INPUT_DIR}" ]]; then
    log_error "Input directory does not exist: ${INPUT_DIR}"
    exit 1
fi

# Validate boolean flags
for flag in SINGLE_END DISABLE_ADAPTER_TRIMMING HTML_REPORT JSON_REPORT OVERWRITE; do
    if ! [[ "${!flag}" =~ ^[tf]$ ]]; then
        log_error "Flag --$(echo "${flag}" | tr 'A-Z_' 'a-z-') must be 't' or 'f' (got '${!flag}')."
        exit 1
    fi
done

# Validate numeric parameters
for num in NSLOTS MIN_LENGTH QUALIFIED_QUALITY_PHRED UNQUALIFIED_PERCENT_LIMIT; do
    if ! [[ "${!num}" =~ ^[0-9]+$ ]]; then
        log_error "Parameter ${num} must be numeric (got '${!num}')."
        exit 1
    fi
done

###############################################################################
# 6. Check dependencies
###############################################################################

log "Checking dependencies..."

check_cmd fastp

###############################################################################
# 7. Prepare output directory
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
    if [[ "${OVERWRITE}" == "t" ]]; then
        log_warn "Overwriting existing directory: ${OUTPUT_DIR}"
        rm -rf "${OUTPUT_DIR}"
    else
        log_error "Output directory exists: ${OUTPUT_DIR}. Use --overwrite=t to overwrite."
        exit 1
    fi
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}/stats"

###############################################################################
# 8. Find input files
###############################################################################

log "Searching for input files..."

# Use arrays to hold file paths
readarray -t R1_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*${R1_PATTERN}" | sort)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    log_error "No files found with pattern: *${R1_PATTERN}"
    exit 1
fi

if [[ "${SINGLE_END}" == "f" ]]; then
    readarray -t R2_FILES < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*${R2_PATTERN}" | sort)

    if [[ ${#R2_FILES[@]} -eq 0 ]]; then
        log_error "No R2 files found with pattern: *${R2_PATTERN}"
        exit 1
    fi

    if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
        log_error "Number of R1 and R2 files do not match (R1: ${#R1_FILES[@]}, R2: ${#R2_FILES[@]})"
        exit 1
    fi

    log "Found ${#R1_FILES[@]} sample pairs"
else
    log "Found ${#R1_FILES[@]} single-end samples"
fi

###############################################################################
# 9. Process samples
###############################################################################

SUMMARY_FILE="${OUTPUT_DIR}/stats/summary.tsv"
echo -e "sample\ttotal_reads_before\ttotal_bases_before\tq20_bases_before\tq30_bases_before\ttotal_reads_after\ttotal_bases_after\tq20_bases_after\tq30_bases_after\tpercent_passed" > "${SUMMARY_FILE}"

for i in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[$i]}"

    # Extract sample name by removing pattern
    SAMPLE_NAME=$(basename "${R1}" "${R1_PATTERN}")

    log "Processing sample ${SAMPLE_NAME} ($(( i + 1 ))/${#R1_FILES[@]})..."

    # Define output files
    HTML_OUT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.html"
    JSON_OUT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.json"

    # Build fastp command (report only mode)
    FASTP_CMD=(
        fastp
        -i "${R1}"
        -w "${NSLOTS}"
        --disable_quality_filtering
        --disable_length_filtering
        --disable_trim_poly_g
    )

    if [[ "${SINGLE_END}" == "f" ]]; then
        R2="${R2_FILES[$i]}"
        FASTP_CMD+=(-I "${R2}")
    fi

    # Add HTML report if requested
    if [[ "${HTML_REPORT}" == "t" ]]; then
        FASTP_CMD+=(--html "${HTML_OUT}")
    else
        FASTP_CMD+=(--html /dev/null)
    fi

    # Add JSON report if requested
    if [[ "${JSON_REPORT}" == "t" ]]; then
        FASTP_CMD+=(--json "${JSON_OUT}")
    else
        FASTP_CMD+=(--json /dev/null)
    fi

    # Disable adapter trimming if requested
    if [[ "${DISABLE_ADAPTER_TRIMMING}" == "t" ]]; then
        FASTP_CMD+=(--disable_adapter_trimming)
    fi

    # Don't output processed reads (report only mode)
    # Create temporary output file(s)
    TEMP_OUT_R1="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_temp_R1.fastq.gz"
    FASTP_CMD+=(-o "${TEMP_OUT_R1}")

    if [[ "${SINGLE_END}" == "f" ]]; then
        TEMP_OUT_R2="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_temp_R2.fastq.gz"
        FASTP_CMD+=(-O "${TEMP_OUT_R2}")
    fi

    # Run fastp
    if ! "${FASTP_CMD[@]}" 2>&1 | tee "${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.log"; then
        log_error "fastp failed for sample ${SAMPLE_NAME}"
        rm -f "${TEMP_OUT_R1}"
        [[ "${SINGLE_END}" == "f" ]] && rm -f "${TEMP_OUT_R2}"
        exit 1
    fi

    # Remove temporary output files after successful completion
    rm -f "${TEMP_OUT_R1}"
    [[ "${SINGLE_END}" == "f" ]] && rm -f "${TEMP_OUT_R2}"
    
    # Extract summary statistics from JSON if available
    if [[ "${JSON_REPORT}" == "t" && -f "${JSON_OUT}" ]]; then
        # Parse JSON using basic tools
        TOTAL_READS_BEFORE=$(grep -Po '"total_reads":\s*\K[0-9]+' "${JSON_OUT}" | head -1)
        TOTAL_BASES_BEFORE=$(grep -Po '"total_bases":\s*\K[0-9]+' "${JSON_OUT}" | head -1)
        Q20_BASES_BEFORE=$(grep -Po '"q20_bases":\s*\K[0-9]+' "${JSON_OUT}" | head -1)
        Q30_BASES_BEFORE=$(grep -Po '"q30_bases":\s*\K[0-9]+' "${JSON_OUT}" | head -1)
        
        TOTAL_READS_AFTER=$(grep -Po '"total_reads":\s*\K[0-9]+' "${JSON_OUT}" | tail -1)
        TOTAL_BASES_AFTER=$(grep -Po '"total_bases":\s*\K[0-9]+' "${JSON_OUT}" | tail -1)
        Q20_BASES_AFTER=$(grep -Po '"q20_bases":\s*\K[0-9]+' "${JSON_OUT}" | tail -1)
        Q30_BASES_AFTER=$(grep -Po '"q30_bases":\s*\K[0-9]+' "${JSON_OUT}" | tail -1)
        
        # Calculate percent passed
        if [[ ${TOTAL_READS_BEFORE} -gt 0 ]]; then
            PERCENT_PASSED=$(awk "BEGIN {printf \"%.2f\", (${TOTAL_READS_AFTER}/${TOTAL_READS_BEFORE})*100}")
        else
            PERCENT_PASSED="0.00"
        fi
        
        # Write to summary file
        echo -e "${SAMPLE_NAME}\t${TOTAL_READS_BEFORE}\t${TOTAL_BASES_BEFORE}\t${Q20_BASES_BEFORE}\t${Q30_BASES_BEFORE}\t${TOTAL_READS_AFTER}\t${TOTAL_BASES_AFTER}\t${Q20_BASES_AFTER}\t${Q30_BASES_AFTER}\t${PERCENT_PASSED}" >> "${SUMMARY_FILE}"
    fi
    
    log "Completed processing ${SAMPLE_NAME}"
done

###############################################################################
# 10. Generate summary report
###############################################################################

log "Generating summary report..."

SUMMARY_TXT="${OUTPUT_DIR}/summary_report.txt"

cat > "${SUMMARY_TXT}" << EOF
================================================================================
Quality Check Report - fastp
================================================================================
Date: $(date)
Input directory: ${INPUT_DIR}
Output directory: ${OUTPUT_DIR}
Number of samples: ${#R1_FILES[@]}

Parameters:
-----------
Read type: $([ "${SINGLE_END}" == "t" ] && echo "single-end" || echo "paired-end")
Pattern R1: ${R1_PATTERN}
$([ "${SINGLE_END}" == "f" ] && echo "Pattern R2: ${R2_PATTERN}")
Threads: ${NSLOTS}
Mode: Report only (no filtering applied)
Adapter trimming: $([ "${DISABLE_ADAPTER_TRIMMING}" == "f" ] && echo "enabled" || echo "disabled")

Output files:
-------------
- Individual HTML reports: ${OUTPUT_DIR}/reports/*_fastp.html
- Individual JSON reports: ${OUTPUT_DIR}/reports/*_fastp.json
- Processing logs: ${OUTPUT_DIR}/reports/*_fastp.log
- Summary statistics: ${OUTPUT_DIR}/stats/summary.tsv
- This report: ${OUTPUT_DIR}/summary_report.txt

================================================================================
EOF

if [[ -f "${SUMMARY_FILE}" ]]; then
    echo -e "\nSummary Statistics:\n" >> "${SUMMARY_TXT}"
    column -t -s $'\t' "${SUMMARY_FILE}" >> "${SUMMARY_TXT}"
fi

cat "${SUMMARY_TXT}"

###############################################################################
# 11. End
###############################################################################

log "${GREEN}quality_check_fastp.sh completed successfully${NC}"
log "Reports available in: ${OUTPUT_DIR}/reports/"
log "Summary statistics: ${SUMMARY_FILE}"
