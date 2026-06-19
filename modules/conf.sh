###############################################################################
# conf.sh
# Configuration file for MG-Proc pipeline
###############################################################################


###############################################################################
# Load conda environment
###############################################################################

# Check if MG-Proc environment is activated, if not, activate it
if [[ "${CONDA_DEFAULT_ENV}" != "MG-Proc" ]]; then
  echo "Activating MG-Proc conda environment..."
  
  # Try to find conda/mamba
  if command -v mamba &> /dev/null; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate MG-Proc
  elif command -v conda &> /dev/null; then
    eval "$(conda shell hook --shell bash)"
    conda activate MG-Proc
  else
    echo "ERROR: Neither mamba nor conda found in PATH"
    echo "Please install mamba or conda and ensure it's in your PATH"
    exit 1
  fi
  
  # Verify activation was successful
  if [[ "${CONDA_DEFAULT_ENV}" != "MG-Proc" ]]; then
    echo "ERROR: Failed to activate MG-Proc environment"
    echo "Please ensure the environment exists. Create it with:"
    echo "  mamba env create -f environment.yml"
    exit 1
  fi
  
  echo "Environment activated successfully"
fi

###############################################################################
# dirs
###############################################################################

# Get the directory where this script is located
MODULES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# tools (using conda environment)
bbduk="bbduk.sh"
bbmerge="bbmerge.sh"
reformat="reformat.sh"
repair_sh="repair.sh"
seqtk="seqtk"
pear="pear"
bzip2="bzip2"
gunzip="gunzip"
rscript="Rscript"
pigz="pigz"
fq2fa="${MODULES_DIR}/resources/fq2fa.sh"
plots="${MODULES_DIR}/resources/plots.R"
quality_check="${MODULES_DIR}/1.2-quality_check.R"
megahit="megahit"
bwa="bwa"
samtools="samtools"
picard="picard" # Assuming picard is available as a command (e.g., via conda), otherwise this should point to the jar file

###############################################################################
# files
###############################################################################

# BBMap adapters file location in conda environment
ADAPTERS=$(find "${CONDA_PREFIX}/opt/bbmap"* -name "adapters.fa")
if [[ ! -f "${ADAPTERS}" ]]; then
  echo "ERROR: BBMap adapters file not found at expected location: ${ADAPTERS}"
  echo "Please ensure BBMap is installed in the MG-Proc mamba environment."
fi

###############################################################################
# functions
###############################################################################

# Color codes for logging

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

log()        { echo -e "[INFO] $*"; }
log_warn()   { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
log_error()  { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Function to check for required tools and dependencies

function check_dependencies {
  local missing_tools=()
  
  # Check command-line tools
  for tool in bbduk.sh bbmerge.sh reformat.sh repair.sh seqtk pear bzip2 gunzip pigz megahit bwa samtools; do
    if ! command -v ${tool} &> /dev/null; then
      missing_tools+=("${tool}")
    fi
  done
  
  # Check for picard (Java jar file)
  if ! command -v picard &> /dev/null; then
    missing_tools+=("picard")
  fi
  
  # Check for emboss tools (infoseq is part of emboss)
  if ! command -v infoseq &> /dev/null; then
    missing_tools+=("emboss")
  fi
  
  # Check R packages
  if command -v Rscript &> /dev/null; then
    # Check tidyverse
    if ! Rscript -e "library(tidyverse)" &> /dev/null; then
      missing_tools+=("r-tidyverse")
    fi
    # Check ShortRead (Bioconductor)
    if ! Rscript -e "library(ShortRead)" &> /dev/null; then
      missing_tools+=("bioconductor-shortread")
    fi
    # Check doParallel
    if ! Rscript -e "library(doParallel)" &> /dev/null; then
      missing_tools+=("r-doparallel")
    fi
    # Check dada2 (Bioconductor)
    if ! Rscript -e "library(dada2)" &> /dev/null; then
      missing_tools+=("bioconductor-dada2")
    fi
    # Check optparse
    if ! Rscript -e "library(optparse)" &> /dev/null; then
      missing_tools+=("r-optparse")
    fi
  else
    missing_tools+=("R")
  fi
  
  # Check local resources
  if [[ ! -f "${fq2fa}" ]]; then
    missing_tools+=("fq2fa.sh")
  fi
  if [[ ! -f "${plots}" ]]; then
    missing_tools+=("plots.R")
  fi
  if [[ ! -f "${quality_check}" ]]; then
    missing_tools+=("1.2-quality_check.R")
  fi
  
  # Report missing tools
  if [[ ${#missing_tools[@]} -gt 0 ]]; then
    log_error "The following required tools are not installed or not in PATH:"
    printf '  - %s\n' "${missing_tools[@]}" >&2
    log_error ""
    log_error "Please activate the MG-Proc mamba environment:"
    log_error "  mamba activate MG-Proc"
    return 1
  fi
  
  return 0
}

# Function to check that a command exists in PATH
check_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log_error "Required command '$1' not found in PATH."
    exit 1
  fi
}

# Function to check that a file exists
check_file() {
  local f="$1"
  local label="${2:-file}"
  if [[ ! -f "$f" ]]; then
    log_error "$label '$f' does not exist or is not a regular file."
    exit 1
  fi
}


# Function to count the number of reads in a FASTQ file
function count_fastq {

  N=$(wc -l  "${1}" | cut -f1 -d" ")
  N=$( echo "${N}" / 4 | bc -l )
  echo "${N}"

}

# Function to count the number of sequences in a FASTA file
function count_fasta {

  N=$(egrep -c ">" "${1}")
  echo "${N}"

}

# Function to compute the mean length of sequences in a FASTA or FASTQ file
function compute_mean_length {

  N=$(infoseq "${1}" | \
      awk 'NR > 1 {tot = $6 + tot; n++} END {print tot/n}')
  echo "${N}"

}

