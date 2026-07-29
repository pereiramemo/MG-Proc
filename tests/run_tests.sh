#!/usr/bin/env bash
###############################################################################
# run_tests.sh — script-level smoke tests for the MG-Proc bin/ scripts.
#
# Exercises the ported Python scripts directly (no containers), so the pipeline
# tools must be on PATH — e.g. run inside a conda/mamba env that provides fastp,
# bbmap, seqtk, pear, pigz, megahit, bwa, samtools, picard, and matplotlib.
# For the containerized end-to-end run, use `nextflow run mg-proc.nf` instead.
###############################################################################

set -uo pipefail

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${WORKDIR}/bin"
DATA="${WORKDIR}/tests/data"
OUT="${WORKDIR}/tests/outputs"

SE_R1="${DATA}/SRR4831661/SRR4831661.fastq"
PE_R1="${DATA}/SRR12479690/SRR12479690_1.fastq.gz"
PE_R2="${DATA}/SRR12479690/SRR12479690_2.fastq.gz"

pass=0
fail=0

run_test() {
    local name="$1"; shift
    echo "------------------------------------------------------------"
    echo "TEST: ${name}"
    echo "------------------------------------------------------------"
    if "$@"; then
        echo "[PASS] ${name}"; (( pass++ )) || true
    else
        echo "[FAIL] ${name}"; (( fail++ )) || true
    fi
    echo
}

###############################################################################
# 0. Static syntax checks
###############################################################################

run_test "pycompile" python -m py_compile \
    "${BIN}/utils.py" "${BIN}/1.1-quality-check.py" \
    "${BIN}/2-preprocess.py" "${BIN}/3-assembly-and-map.py"

###############################################################################
# 1. Single-end preprocessing + assembly
###############################################################################

run_test "se_preprocess" \
    "${BIN}/2-preprocess.py" \
        --reads "${SE_R1}" \
        --single_end t \
        --reformat t \
        --subsample t \
        --plot t \
        --sample_name se_test \
        --output_dir "${OUT}/se_preprocess" \
        --overwrite t

run_test "se_assembly" \
    "${BIN}/3-assembly-and-map.py" \
        --reads1 "${OUT}/se_preprocess/output/se_test_se_qc-02.fastq" \
        --single_end t \
        --sample_name se_test \
        --output_dir "${OUT}/se_assembly" \
        --overwrite t

###############################################################################
# 2. Paired-end preprocessing + assembly
###############################################################################

run_test "pe_preprocess" \
    "${BIN}/2-preprocess.py" \
        --reads "${PE_R1}" \
        --reads2 "${PE_R2}" \
        --subsample t \
        --repair t \
        --output_pe t \
        --plot t \
        --sample_name pe_test \
        --output_dir "${OUT}/pe_preprocess" \
        --overwrite t

run_test "pe_assembly" \
    "${BIN}/3-assembly-and-map.py" \
        --reads1 "${OUT}/pe_preprocess/output/pe_test_R1_qc-02.fastq" \
        --reads2 "${OUT}/pe_preprocess/output/pe_test_R2_qc-02.fastq" \
        --sample_name pe_test \
        --output_dir "${OUT}/pe_assembly" \
        --overwrite t

###############################################################################
# 3. fastp QC report (single-end and paired-end)
###############################################################################

run_test "se_quality_check" \
    "${BIN}/1.1-quality-check.py" \
        --reads1 "${SE_R1}" \
        --single_end t \
        --sample_name se_test \
        --output_dir "${OUT}/se_quality_check" \
        --overwrite t

run_test "pe_quality_check" \
    "${BIN}/1.1-quality-check.py" \
        --reads1 "${PE_R1}" \
        --reads2 "${PE_R2}" \
        --sample_name pe_test \
        --output_dir "${OUT}/pe_quality_check" \
        --overwrite t

###############################################################################
# Summary
###############################################################################

echo "============================================================"
echo "Results: ${pass} passed, ${fail} failed"
echo "============================================================"

[[ ${fail} -eq 0 ]]
