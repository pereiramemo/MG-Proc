#!/usr/bin/env bash
###############################################################################
### 3-assembly_and_map_pipeline.sh
### De novo assembly + read mapping + duplicate removal
###############################################################################

set -euo pipefail

###############################################################################
### 1. Set env
###############################################################################

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load configuration (must define paths to tools: megahit, bwa, samtools, picard, etc.)
source "${SCRIPT_DIR}/conf.sh"

###############################################################################
### 1.1 Helper function
###############################################################################

# Find contigs file in assembly directory
# Searches for contigs files with various extensions in multiple locations
# Supports both compressed (.gz) and uncompressed files

find_contigs_file() {
  local assem_dir="$1"
  local sample_name="$2"
  local extensions=("contigs.fa" "contigs.fasta" "contigs.fna" "fa" "fasta" "fna" \
                    "contigs.fa.gz" "contigs.fasta.gz" "contigs.fna.gz" "fa.gz" "fasta.gz" "fna.gz")
  local search_paths=(
    "${assem_dir}/${sample_name}"
    "${assem_dir}"
  )

  # First, look for files with sample name prefix
  for dir in "${search_paths[@]}"; do
    for ext in "${extensions[@]}"; do
      local candidate="${dir}/${sample_name}.${ext}"
      if [[ -f "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    done
  done

  # If not found, look for any .contigs.{fa,fasta,fna}[.gz] file in the directories
  # Enable nullglob so unmatched patterns expand to nothing
  shopt -s nullglob
  for dir in "${search_paths[@]}"; do
    for ext in fa fasta fna; do
      # Check uncompressed files
      for candidate in "${dir}"/*.contigs."${ext}"; do
        if [[ -f "${candidate}" ]]; then
          log_warn "Found contigs file '${candidate}' but it doesn't match expected pattern '${sample_name}.contigs.*'"
          log_warn "Using it anyway. Consider renaming or using --contigs to specify explicitly."
          shopt -u nullglob  # Disable nullglob before returning
          echo "${candidate}"
          return 0
        fi
      done
      # Check compressed files
      for candidate in "${dir}"/*.contigs."${ext}".gz; do
        if [[ -f "${candidate}" ]]; then
          log_warn "Found contigs file '${candidate}' but it doesn't match expected pattern '${sample_name}.contigs.*'"
          log_warn "Using it anyway. Consider renaming or using --contigs to specify explicitly."
          shopt -u nullglob  # Disable nullglob before returning
          echo "${candidate}"
          return 0
        fi
      done
    done
  done
  shopt -u nullglob  # Disable nullglob before returning

  # Not found
  return 1
}

###############################################################################
### 2. Define help
###############################################################################

show_usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --reads1 CHAR              Input R1 (or single-end) metagenome reads (fastq/fa)
  --sample_name CHAR         Sample name used to name output files

Paired-end only (ignored when --single_end=t):
  --reads2 CHAR              Input R2 metagenome reads (fastq/fa)

Optional:
  --single_end t|f           Process as single-end reads (default: f)
  --contigs CHAR             Path to pre-assembled contigs file (FASTA format)
                             Supports both compressed (.gz) and uncompressed files
                             Takes precedence over --assem_dir
  --assem_dir CHAR           Directory with previously computed assemblies
                             Will search for: SAMPLE_NAME.contigs.{fa,fasta,fna}[.gz]
                             in ASSEM_DIR/ or ASSEM_DIR/SAMPLE_NAME/
                             Supports both compressed (.gz) and uncompressed files
  --assem_preset CHAR        MEGAHIT preset to generate assembly
                             (default: meta-sensitive)
  --nslots NUM               Number of threads used (default: 12)
  --min_contig_length NUM    Minimum length of contigs to keep (default: 250)
  --output_dir CHAR          Output directory (default: mg-clust_output-1)
  --overwrite t|f            Overwrite previous folder if present (default: f)
  --remove_duplicates t|f    Remove PCR duplicates with Picard (default: f)
  --help                     Print this help and exit
EOF
}

###############################################################################
### 3. Parse input parameters (using getopt)
###############################################################################

# Defaults (can be overridden by CLI)
CONTIGS=""
ASSEM_DIR=""
ASSEM_PRESET="meta-sensitive"
NSLOTS=12
MIN_CONTIG_LENGTH=250
OUTPUT_DIR="mg-clust_output-1"
OVERWRITE="f"
REMOVE_DUPLICATES="f"
SINGLE_END="f"
R1=""
R2=""
SAMPLE_NAME=""

# Ensure getopt is available
check_cmd getopt

ARGS=$(
  getopt -o '' \
    --long help,contigs:,assem_dir:,assem_preset:,nslots:,min_contig_length:,output_dir:,overwrite:,remove_duplicates:,single_end:,reads1:,reads2:,sample_name: \
    -n "$(basename "$0")" -- "$@" \
) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    --help)
      show_usage
      exit 0
      ;;
    --contigs)
      CONTIGS="$2"; shift 2 ;;
    --assem_dir)
      ASSEM_DIR="$2"; shift 2 ;;
    --assem_preset)
      ASSEM_PRESET="$2"; shift 2 ;;
    --nslots)
      NSLOTS="$2"; shift 2 ;;
    --min_contig_length)
      MIN_CONTIG_LENGTH="$2"; shift 2 ;;
    --output_dir)
      OUTPUT_DIR="$2"; shift 2 ;;
    --overwrite)
      OVERWRITE="$2"; shift 2 ;;
    --remove_duplicates)
      REMOVE_DUPLICATES="$2"; shift 2 ;;
    --single_end)
      SINGLE_END="$2"; shift 2 ;;
    --reads1)
      R1="$2"; shift 2 ;;
    --reads2)
      R2="$2"; shift 2 ;;
    --sample_name)
      SAMPLE_NAME="$2"; shift 2 ;;
    --)
      shift
      break
      ;;
    *)
      log_error "Unexpected argument: $1"
      show_usage
      exit 1
      ;;
  esac
done

###############################################################################
### 4 Validate parameters
###############################################################################

# Mandatory parameters
for v in R1 SAMPLE_NAME; do
  if [[ -z "${!v}" ]]; then
    log_error "Missing mandatory parameter: --$(echo "$v" | tr '[:upper:]' '[:lower:]')"
    show_usage
    exit 1
  fi
done

if [[ "${SINGLE_END}" == "f" && -z "${R2}" ]]; then
  log_error "Missing mandatory parameter: --reads2 (required for paired-end mode)"
  show_usage
  exit 1
fi

# Simple validation for boolean flags
for flag in OVERWRITE REMOVE_DUPLICATES SINGLE_END; do
  if ! [[ "${!flag}" =~ ^[tf]$ ]]; then
    log_error "--$(echo "${flag}" | tr '[:upper:]' '[:lower:]') must be 't' or 'f' (got '${!flag}')"
    exit 1
  fi
done

# Numeric validations
if ! [[ "${NSLOTS}" =~ ^[0-9]+$ ]]; then
  log_error "--nslots must be an integer (got '${NSLOTS}')"
  exit 1
fi

if ! [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]]; then
  log_error "--min_contig_length must be an integer (got '${MIN_CONTIG_LENGTH}')"
  exit 1
fi

###############################################################################
### 5. Check tools
###############################################################################

check_dependencies

###############################################################################
### 6. Check mandatory files
###############################################################################

check_file "${R1}" "read1"
if [[ "${SINGLE_END}" == "f" ]]; then
  check_file "${R2}" "read2"
fi

###############################################################################
### 7. Create output folder
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${OVERWRITE}" == "t" ]]; then
    log_warn "Output directory '${OUTPUT_DIR}' exists and will be overwritten."
    rm -rf "${OUTPUT_DIR}" || {
      log_error "Failed to remove existing output directory '${OUTPUT_DIR}'."
      exit 1
    }
  else
    log_error "Output directory '${OUTPUT_DIR}' already exists; use --overwrite t to overwrite."
    exit 1
  fi
fi

# Create base output directory if assembly is not run by megahit (when CONTIGS or ASSEM_DIR is given)
if [[ -n "${CONTIGS}" || -n "${ASSEM_DIR}" ]]; then
  mkdir -p "${OUTPUT_DIR}" || {
    log_error "Failed to create output directory '${OUTPUT_DIR}'."
    exit 1
  }
fi

# Set up trap to clean up tmp directory on exit (success or failure)
trap 'rm -rf "${OUTPUT_DIR}/tmp" 2>/dev/null' EXIT

###############################################################################
### 8. De novo assembly or identify external contigs
###############################################################################

ASSEMBLY_FILE=""
ASSEMBLY_FILE_LOCAL=""  # Track if we created a local copy of an external assembly

# Priority 1: Use explicitly provided contigs file
if [[ -n "${CONTIGS}" ]]; then
  log "Using externally provided contigs file: ${CONTIGS}"
  ASSEMBLY_FILE="${CONTIGS}"
  check_file "${ASSEMBLY_FILE}" "contigs file"

# Priority 2: Search for contigs in provided assembly directory
elif [[ -n "${ASSEM_DIR}" ]]; then
  log "Searching for contigs file in assembly directory: ${ASSEM_DIR}"

  # Validate assembly directory exists
  if [[ ! -d "${ASSEM_DIR}" ]]; then
    log_error "Assembly directory '${ASSEM_DIR}' does not exist or is not a directory."
    exit 1
  fi

  # Try to find contigs file
  ASSEMBLY_FILE=$(find_contigs_file "${ASSEM_DIR}" "${SAMPLE_NAME}")

  if [[ -z "${ASSEMBLY_FILE}" ]]; then
    log_error "Could not find contigs file for sample '${SAMPLE_NAME}' in '${ASSEM_DIR}'."
    log_error "Searched in:"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.{contigs.fa,contigs.fasta,contigs.fna,fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}.{contigs.fa,contigs.fasta,contigs.fna,fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}/*.contigs.{fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/*.contigs.{fa,fasta,fna}[.gz]"
    log_error ""
    log_error "Use --contigs to specify the exact path to your contigs file."
    exit 1
  fi

  log "Found contigs file: ${ASSEMBLY_FILE}"

# Priority 3: Run de novo assembly with MEGAHIT
else
  log "Running MEGAHIT assembly..."

  if [[ "${SINGLE_END}" == "t" ]]; then
    "${megahit}" \
      --num-cpu-threads "${NSLOTS}" \
      -r "${R1}" \
      --presets "${ASSEM_PRESET}" \
      --min-contig-len "${MIN_CONTIG_LENGTH}" \
      --out-prefix "${SAMPLE_NAME}" \
      --out-dir "${OUTPUT_DIR}"
  else
    "${megahit}" \
      --num-cpu-threads "${NSLOTS}" \
      -1 "${R1}" \
      -2 "${R2}" \
      --presets "${ASSEM_PRESET}" \
      --min-contig-len "${MIN_CONTIG_LENGTH}" \
      --out-prefix "${SAMPLE_NAME}" \
      --out-dir "${OUTPUT_DIR}"
  fi

  log "MEGAHIT completed."

  ASSEMBLY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.contigs.fa"
  check_file "${ASSEMBLY_FILE}" "assembly file"
fi

# Resolve symbolic links to avoid issues with BWA index file locations
if [[ -L "${ASSEMBLY_FILE}" ]]; then
  ASSEMBLY_FILE_ORIGINAL="${ASSEMBLY_FILE}"
  ASSEMBLY_FILE=$(readlink -f "${ASSEMBLY_FILE}")
  log "Resolved symbolic link: ${ASSEMBLY_FILE_ORIGINAL} -> ${ASSEMBLY_FILE}"
fi

# Handle external assembly files (copy to OUTPUT_DIR to avoid polluting original directory with BWA indices)
if [[ -n "${CONTIGS}" || -n "${ASSEM_DIR}" ]]; then

  # Check if file is compressed
  if [[ "${ASSEMBLY_FILE}" == *.gz ]]; then
    log "Detected gzipped assembly file. Decompressing to output directory..."
    ASSEMBLY_FILE_LOCAL="${OUTPUT_DIR}/$(basename "${ASSEMBLY_FILE}" .gz)"

    "${gunzip}" -c "${ASSEMBLY_FILE}" > "${ASSEMBLY_FILE_LOCAL}" || {
      log_error "Failed to decompress assembly file."
      exit 1
    }

    log "Decompressed: $(basename "${ASSEMBLY_FILE}") -> $(basename "${ASSEMBLY_FILE_LOCAL}")"
  else
    log "Copying external assembly file to output directory..."
    ASSEMBLY_FILE_LOCAL="${OUTPUT_DIR}/$(basename "${ASSEMBLY_FILE}")"

    cp "${ASSEMBLY_FILE}" "${ASSEMBLY_FILE_LOCAL}" || {
      log_error "Failed to copy assembly file to output directory."
      exit 1
    }

    log "Copied: $(basename "${ASSEMBLY_FILE}") -> output directory"
  fi

  # Point ASSEMBLY_FILE to the local copy for downstream processing
  ASSEMBLY_FILE="${ASSEMBLY_FILE_LOCAL}"
  check_file "${ASSEMBLY_FILE}" "assembly file in output directory"
fi

###############################################################################
### 9. Map short reads
###############################################################################

log "Checking number of assembled contigs in '${ASSEMBLY_FILE}'..."
ASSEMBLY_FILE_NSEQ=$(grep -c ">" "${ASSEMBLY_FILE}" || echo 0)

if (( ASSEMBLY_FILE_NSEQ < 5 )); then 
  log_warn "Not enough assembled sequences to continue (${ASSEMBLY_FILE_NSEQ} < 5). Exiting gracefully."
  # Distinguish from success but not a hard error
  exit 0
fi

log "Found ${ASSEMBLY_FILE_NSEQ} contigs. Proceeding with mapping."

# Index assembly
log "Indexing assembly with BWA..."
"${bwa}" index "${ASSEMBLY_FILE}" || {
  log_error "bwa index failed."
  exit 1
}

# Map reads and convert directly to BAM (avoid huge SAM files)
log "Mapping reads with BWA-MEM and converting to BAM (q>=10, primary alignments)..."
if [[ "${SINGLE_END}" == "t" ]]; then
  "${bwa}" mem -M -t "${NSLOTS}" "${ASSEMBLY_FILE}" "${R1}" | \
    "${samtools}" view -@ "${NSLOTS}" -q 10 -F 260 -b > "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
    log_error "bwa mem or samtools view failed."
    exit 1
  }
else
  "${bwa}" mem -M -t "${NSLOTS}" "${ASSEMBLY_FILE}" "${R1}" "${R2}" | \
    "${samtools}" view -@ "${NSLOTS}" -q 10 -F 260 -b > "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
    log_error "bwa mem or samtools view failed."
    exit 1
  }
fi

# Sort BAM
log "Sorting BAM..."
"${samtools}" sort \
  -@ "${NSLOTS}" \
  -o "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
  "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "samtools sort failed."
  exit 1
}

# Index sorted BAM
log "Indexing sorted BAM..."
"${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" || {
  log_error "samtools index on sorted BAM failed."
  exit 1
}

# Remove duplicates with Picard (optional)
if [[ "${REMOVE_DUPLICATES}" == "t" ]]; then
  log "Marking and removing duplicates with Picard..."
  mkdir -p "${OUTPUT_DIR}/tmp" || {
    log_error "Failed to create temporary directory '${OUTPUT_DIR}/tmp'."
    exit 1
  }

  java -jar "${picard}" MarkDuplicates \
    INPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
    OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam" \
    METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.metrics.txt" \
    REMOVE_DUPLICATES=TRUE \
    ASSUME_SORTED=TRUE \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
    TMP_DIR="${OUTPUT_DIR}/tmp" || {
    log_error "Picard MarkDuplicates failed."
    exit 1
  }

  # Index final markdup BAM
  log "Indexing final duplicate-removed BAM..."
  "${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam" || {
    log_error "samtools index on markdup BAM failed."
    exit 1
  }

  FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam"
else
  log "Skipping duplicate removal (--remove_duplicates=f)"
  FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"
fi

###############################################################################
### 10. Clean
###############################################################################

log "Removing intermediate mapping files..."
# Always remove the initial unsorted BAM
rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "Failed to remove intermediate BAM file."
  exit 1
}

# If duplicates were removed, also remove the sorted BAM (it's intermediate)
if [[ "${REMOVE_DUPLICATES}" == "t" ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
        "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam.bai" || {
    log_warn "Failed to remove intermediate sorted BAM files (non-fatal)."
  }
fi

# Remove BWA index files (always cleanup since they were created by this script)
log "Removing BWA index files..."
rm -f "${ASSEMBLY_FILE}".{amb,ann,bwt,pac,sa} || {
  log_warn "Failed to remove some BWA index files (non-fatal)."
}

# Remove temporary local copy of external assembly file if created
if [[ -n "${ASSEMBLY_FILE_LOCAL}" && -f "${ASSEMBLY_FILE_LOCAL}" ]]; then
  log "Removing temporary assembly file from output directory..."
  rm -f "${ASSEMBLY_FILE_LOCAL}" || {
    log_warn "Failed to remove temporary assembly file (non-fatal)."
  }
fi

# Note: tmp directory is cleaned up automatically by trap on exit

if [[ -z "${CONTIGS}" && -z "${ASSEM_DIR}" ]]; then
  # Remove MEGAHIT intermediate contigs but keep main assembly
  if [[ -d "${OUTPUT_DIR}/intermediate_contigs" ]]; then
    log "Removing MEGAHIT intermediate contigs..."
    rm -rf "${OUTPUT_DIR}/intermediate_contigs" || {
      log_warn "Failed to remove 'intermediate_contigs' directory (non-fatal)."
    }
  fi
fi

###############################################################################
### 11. Exit
###############################################################################

log "${GREEN}assembly_and_map_pipeline.sh exited successfully.${NC}"
log "Final BAM file: ${FINAL_BAM}"
exit 0
