#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Import shared helpers from bin/utils.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import log, log_warn, log_error, derive_sample_name, build_log, decompress_or_link

SCRIPT_NAME = "1.1-quality-check.py"
SCRIPT_DESC = "Run a fastp QC report on a single sample (report only; no filtering applied)."

# Dev only section for testing in IDE (uncomment to test)
""" 
reads1 = "/home/epereira/workspace/repos/tools/MG-Proc/tests/data/SRR12479690/SRR12479690_1.fastq.gz"
reads2 = "/home/epereira/workspace/repos/tools/MG-Proc/tests/data/SRR12479690/SRR12479690_2.fastq.gz"
output_dir = Path("/home/epereira/workspace/repos/tools/MG-Proc/tests/output/1.1-quality-check-out")
sample_name = "SRR12479690"
single_end = False
nslots = 12
min_length = 50
qualified_quality_phred = 20
unqualified_percent_limit = 40
disable_adapter_trimming = True
html_report = True
json_report = True
overwrite = False
"""

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Run a fastp QC report on a single sample (report only; no filtering applied)"
    )
    p.add_argument("--reads1",       required=True,  help="R1 FASTQ file, or single-end file when --single_end t (required)")
    p.add_argument("--reads2",       default=None,   help="R2 FASTQ file (paired-end only)")
    p.add_argument("--output_dir",   required=True,  help="Output directory (required)")
    p.add_argument("--sample_name",  default=None,   help="Sample name [default: derived from R1 filename]")
    p.add_argument("--single_end",   choices=["t", "f"], default="f", help="Process single-end reads [default=f]")
    p.add_argument("--nslots",       type=int,   default=12,  help="Number of threads [default=12]")
    p.add_argument("--min_length",   type=int,   default=50,  help="Minimum read length filter, reporting only [default=50]")
    p.add_argument("--qualified_quality_phred",   type=int, default=20, help="Minimum quality value for a qualified base, reporting only [default=20]")
    p.add_argument("--unqualified_percent_limit", type=int, default=40, help="Maximum percent of unqualified bases allowed, reporting only [default=40]")
    p.add_argument("--disable_adapter_trimming",  choices=["t", "f"], default="t", help="Disable adapter trimming in report [default=t]")
    p.add_argument("--html_report",  choices=["t", "f"], default="t", help="Generate HTML report [default=t]")
    p.add_argument("--json_report",  choices=["t", "f"], default="t", help="Keep JSON report after stats extraction [default=t]")
    p.add_argument("--overwrite",    choices=["t", "f"], default="f", help="Overwrite existing output directory [default=f]")
    return p.parse_args()

################################################################################
# 3. Define main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    reads1                    = opts.reads1
    reads2                    = opts.reads2
    output_dir                = Path(opts.output_dir)
    sample_name               = opts.sample_name
    single_end                = opts.single_end == "t"
    nslots                    = opts.nslots
    min_length                = opts.min_length
    qualified_quality_phred   = opts.qualified_quality_phred
    unqualified_percent_limit = opts.unqualified_percent_limit
    disable_adapter_trimming  = opts.disable_adapter_trimming == "t"
    html_report               = opts.html_report == "t"
    json_report               = opts.json_report == "t"
    overwrite                 = opts.overwrite == "t"

    ###########################################################################
    # Step 2: Validate inputs and tools, create output directories
    ###########################################################################

    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not single_end:
        if reads2 is None:
            log_error("--reads2 is required for paired-end mode.")
            sys.exit(1)
        if not os.path.isfile(reads2):
            log_error(f"R2 file does not exist: {reads2}")
            sys.exit(1)

    log("Checking dependencies...")
    if not shutil.which("fastp"):
        log_error("fastp not found. Please install it or add it to PATH.")
        sys.exit(1)

    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    results_dir = output_dir / "output"
    stats_dir   = output_dir / "stats"
    logs_dir    = output_dir / "logs"
    for d in (results_dir, stats_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    if sample_name is None:
        sample_name = derive_sample_name(reads1)
    log(f"Processing sample: {sample_name}")

    html_out  = results_dir / f"{sample_name}_fastp.html"
    json_out  = results_dir / f"{sample_name}_fastp.json"
    log_out   = logs_dir    / f"1.1-quality-check-{sample_name}.log"
    stats_out = stats_dir   / f"1.1-quality-check-{sample_name}-stats.tsv"

    ###########################################################################
    # Step 3: Build and execute fastp command (report-only mode)
    ###########################################################################

    mode = "single-end" if single_end else "paired-end"
    params = [
        f"Threads: {nslots}",
        f"Read type: {mode}",
        "Mode: report only (no filtering applied)",
        f"Minimum length (reporting): {min_length}",
        f"Qualified quality phred (reporting): {qualified_quality_phred}",
        f"Unqualified percent limit (reporting): {unqualified_percent_limit}",
        f"Adapter trimming: {'disabled' if disable_adapter_trimming else 'enabled'}",
    ]
    inputs = [f"R1: {reads1}"] + ([] if single_end else [f"R2: {reads2}"])

    # fastp only auto-detects gzip; bzip2 (or any other) input must be
    # decompressed first, or it silently misparses the raw bytes as FASTQ.
    with tempfile.TemporaryDirectory(prefix=f"1.1-quality-check-{sample_name}-",
                                     dir=str(output_dir)) as tmp_dir:
        reads1_run = os.path.join(tmp_dir, "R1.fastq")
        decompress_or_link(reads1, reads1_run)
        reads2_run = None
        if not single_end:
            reads2_run = os.path.join(tmp_dir, "R2.fastq")
            decompress_or_link(reads2, reads2_run)

        cmd = [
            "fastp",
            "-i", reads1_run,
            "-w", str(nslots),
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--disable_trim_poly_g",
            "--length_required", str(min_length),
            "--qualified_quality_phred", str(qualified_quality_phred),
            "--unqualified_percent_limit", str(unqualified_percent_limit),
            "--json", str(json_out),
        ]
        cmd += ["-I", reads2_run] if not single_end else []
        cmd += ["--html", str(html_out)] if html_report else ["--html", "/dev/null"]
        cmd += ["--disable_adapter_trimming"] if disable_adapter_trimming else []

        log("Running fastp...")
        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        print(result.stdout, end="")

    ###########################################################################
    # Step 4: Check fastp result
    ###########################################################################

    if result.returncode != 0:
        log_error(f"fastp failed for sample {sample_name}")
        log_out.write_text(build_log(
            script_name=SCRIPT_NAME,
            script_desc=SCRIPT_DESC,
            sample_name=sample_name,
            inputs=inputs, params=params,
            outputs=[f"HTML report: {html_out}", f"JSON report: {json_out}"],
            command=" ".join([SCRIPT_NAME] + sys.argv[1:]),
            exit_status=result.returncode,
            tool_log=result.stdout,
        ))
        sys.exit(1)

    ###########################################################################
    # Step 5: Create stats and log files
    ###########################################################################

    with open(json_out) as f:
        bf = json.load(f)["summary"]["before_filtering"]

    total_reads    = bf["total_reads"]
    total_bases    = bf["total_bases"]
    q20_bases      = bf["q20_bases"]
    q20_rate       = bf["q20_rate"]
    q30_bases      = bf["q30_bases"]
    q30_rate       = bf["q30_rate"]
    r1_mean_length = bf.get("read1_mean_length", "NA")
    r2_mean_length = bf.get("read2_mean_length", "NA")
    gc_content     = bf["gc_content"]

    header = "\t".join([
        "sample", "total_reads", "total_bases",
        "q20_bases", "q20_rate", "q30_bases", "q30_rate",
        "read1_mean_length", "read2_mean_length", "gc_content",
    ])
    row = "\t".join(str(v) for v in [
        sample_name, total_reads, total_bases,
        q20_bases, q20_rate, q30_bases, q30_rate,
        r1_mean_length, r2_mean_length, gc_content,
    ])
    stats_out.write_text(header + "\n" + row + "\n")

    if not json_report:
        json_out.unlink(missing_ok=True)

    log("\033[0;32m1.1-quality-check.py completed successfully\033[0m")

    log_out.write_text(build_log(
        script_name=SCRIPT_NAME,
        script_desc=SCRIPT_DESC,
        sample_name=sample_name,
        inputs=inputs, params=params,
        outputs=[
            f"HTML report: {html_out}",
            f"JSON report: {json_out if json_report else 'not kept'}",
            f"Statistics: {stats_out}",
        ],
        command=" ".join([SCRIPT_NAME] + sys.argv[1:]),
        exit_status=0,
        tool_log=result.stdout,
    ))

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
