#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Import shared helpers from bin/utils.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (log, log_warn, log_error, derive_sample_name, build_log,
                   count_fastq, mean_length, decompress_or_link)

SCRIPT_NAME = "2-preprocess.py"
SCRIPT_DESC = ("Preprocess raw Illumina metagenomic reads (paired-end or single-end): "
               "optional reformat/repair, subsample, adapter trim, quality trim, "
               "PE merging, and FASTA conversion.")

# Dev only section
"""
r1 = Path("/home/epereira/workspace/repos/tools/MG-Proc/tests/data/SRR12479690/SRR12479690_1.fastq.gz")               
r2 = Path("/home/epereira/workspace/repos/tools/MG-Proc/tests/data/SRR12479690/SRR12479690_2.fastq.gz")
output_dir = Path("/home/epereira/workspace/repos/tools/MG-Proc/tests/output/2-preprocess-out")
merger = "pear"
min_overlap = 10
pvalue = 0.01
output_pe = "t"
output_merged = "t"
single_end = "f"
sample_name = "SRR12479690"
"""

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(description=SCRIPT_DESC)
    # Required
    p.add_argument("--reads",        required=True, help="R1 file, or single-end file when --single_end t (required)")
    p.add_argument("--output_dir",   required=True, help="Output directory (required)")
    # Paired-end only
    p.add_argument("--reads2",        default=None,   help="R2 file (paired-end only)")
    p.add_argument("--merger",        choices=["pear", "bbmerge"], default="pear", help="PE merger [default=pear]")
    p.add_argument("--min_overlap",   type=int,   default=10,     help="Minimum PE overlap for PEAR [default=10]")
    p.add_argument("--output_pe",     choices=["t", "f"], default="f", help="Output QC'ed paired-end reads [default=f]")
    p.add_argument("--output_merged", choices=["t", "f"], default="t", help="Output merged QC'ed reads [default=t]")
    p.add_argument("--pvalue",        type=float, default=0.01,    help="p-value for PEAR [default=0.01]")
    # Optional
    p.add_argument("--single_end",    choices=["t", "f"], default="f", help="Process single-end reads [default=f]")
    p.add_argument("--reformat",      choices=["t", "f"], default="f", help="Reformat FASTQ with reformat.sh [default=f]")
    p.add_argument("--repair",        choices=["t", "f"], default="f", help="Reformat + repair FASTQ (PE repair) [default=f]")
    p.add_argument("--clean",         choices=["t", "f"], default="f", help="Remove intermediates [default=f]")
    p.add_argument("--compress",      choices=["t", "f"], default="f", help="Compress outputs with pigz [default=f]")
    p.add_argument("--min_length",    type=int,   default=75,     help="Minimum read length after trimming [default=75]")
    p.add_argument("--min_qual",      type=int,   default=20,     help="Quality trim threshold [default=20]")
    p.add_argument("--nslots",        type=int,   default=12,     help="Threads [default=12]")
    p.add_argument("--sample_name",   default=None, help="Name prefix [default: derived from --reads filename]")
    p.add_argument("--seed",          type=int,   default=123,    help="Random seed for subsampling [default=123]")
    p.add_argument("--subsample",     choices=["t", "f"], default="f", help="Subsample to 10k reads [default=f]")
    p.add_argument("--trim_adapters", choices=["t", "f"], default="f", help="Remove adapters [default=f]")
    p.add_argument("--adapters",      default="", help="Adapters FASTA for BBDuk [default: auto-detect from bbmap]")
    p.add_argument("--overwrite",     choices=["t", "f"], default="f", help="Replace existing output dir [default=f]")
    return p.parse_args()


def find_adapters(explicit):
    """Locate BBDuk's bundled adapters.fa; fall back to the 'adapters' shorthand."""
    if explicit and os.path.isfile(explicit):
        return explicit
    candidates = []
    prefixes = {os.environ.get("CONDA_PREFIX", ""), sys.prefix}
    for base in filter(None, prefixes):
        candidates += glob.glob(f"{base}/opt/bbmap*/resources/adapters.fa")
        candidates += glob.glob(f"{base}/share/bbmap*/resources/adapters.fa")
    bb = shutil.which("bbduk.sh")
    if bb:
        rp = os.path.dirname(os.path.realpath(bb))
        candidates += glob.glob(os.path.join(rp, "resources", "adapters.fa"))
        candidates += glob.glob(os.path.join(rp, "..", "opt", "bbmap*", "resources", "adapters.fa"))
    for c in candidates:
        if os.path.isfile(c):
            return c
    return "adapters"  # BBDuk keyword shorthand for the packaged adapters


def parse_order(basename):
    """Extract the pipeline-step number from a '...-NN.fastq[.gz]' intermediate name."""
    m = re.search(r'-(\d+)\.fastq(?:\.gz)?$', basename)
    return int(m.group(1)) if m else 99


################################################################################
# 3. Define main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    r1            = opts.reads
    r2            = opts.reads2
    output_dir    = Path(opts.output_dir)
    single_end    = opts.single_end == "t"
    reformat      = opts.reformat == "t"
    repair        = opts.repair == "t"
    subsample     = opts.subsample == "t"
    trim_adapters = opts.trim_adapters == "t"
    adapters      = opts.adapters
    output_pe     = opts.output_pe == "t"
    output_merged = opts.output_merged == "t"
    merger        = opts.merger
    min_overlap   = opts.min_overlap
    min_length    = opts.min_length
    min_qual      = opts.min_qual
    pvalue        = opts.pvalue
    seed          = opts.seed
    clean         = opts.clean == "t"
    compress      = opts.compress == "t"
    nslots        = opts.nslots
    sample_name   = opts.sample_name
    if sample_name is None:
        sample_name = derive_sample_name(r1, strip_read_suffix=True)
    overwrite     = opts.overwrite == "t"

    tool_log = []          # accumulated third-party tool output for the log file
    intermediates = []     # paths of intermediate FASTQs (candidates for cleanup)

    ###########################################################################
    # Step 2: Validate inputs and tools
    ###########################################################################

    if not os.path.isfile(r1):
        log_error(f"--reads file does not exist: {r1}")
        sys.exit(1)
    if not single_end:
        if not r2:
            log_error("--reads2 is required for paired-end mode.")
            sys.exit(1)
        if not os.path.isfile(r2):
            log_error(f"--reads2 file does not exist: {r2}")
            sys.exit(1)

    log("Checking dependencies...")
    required = ["bbduk.sh", "seqtk"]
    if reformat or repair:
        required.append("reformat.sh")
    if repair and not single_end:
        required.append("repair.sh")
    if not single_end and output_merged:
        required.append(merger)
    if compress:
        required.append("pigz")
    missing = [t for t in required if not shutil.which(t)]
    if missing:
        log_error("Missing required tools: " + ", ".join(missing))
        sys.exit(1)

    ###########################################################################
    # Step 3: Prepare output directories
    ###########################################################################

    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    results_dir = output_dir / "output"
    logs_dir    = output_dir / "logs"
    stats_dir   = output_dir / "stats"
    for d in (results_dir, logs_dir, stats_dir):
        d.mkdir(parents=True, exist_ok=True)

    log_out   = logs_dir  / f"2-preprocess-{sample_name}.log"
    stats_out = stats_dir / f"2-preprocess-{sample_name}-stats.tsv"
    adapters_found  = find_adapters(adapters)

    ###########################################################################
    # Step 4: Log run parameters
    ###########################################################################

    mode = "single-end" if single_end else "paired-end"
    inputs = [f"R1: {r1}"] + ([] if single_end else [f"R2: {r2}"])
    params = [
        f"Sample name: {sample_name}", f"Read type: {mode}",
        f"Reformat: {'yes' if reformat else 'no'}", f"Repair: {'yes' if repair else 'no'}",
        f"Subsample: {'yes' if subsample else 'no'} (seed {seed})",
        f"Trim adapters: {'yes' if trim_adapters else 'no'}",
        f"Output PE: {'yes' if output_pe else 'no'}", f"Output merged: {'yes' if output_merged else 'no'}",
        f"Merger: {merger}", f"Min overlap: {min_overlap}", f"PEAR p-value: {pvalue}",
        f"Min length: {min_length}", f"Min quality: {min_qual}",
        f"Threads: {nslots}", f"Compress: {'yes' if compress else 'no'}",
        f"Clean: {'yes' if clean else 'no'}",
    ]
    command = " ".join([SCRIPT_NAME] + sys.argv[1:])

    ###########################################################################
    # Step 5: Define helper functions for subprocess execution and logging
    ###########################################################################

    def fail(msg):
        log_error(msg)
        log_out.write_text(build_log(
            script_name=SCRIPT_NAME,
            script_desc=SCRIPT_DESC,
            sample_name=sample_name,
            inputs=inputs,
            params=params,
            outputs=["(failed before completion)"], command=command,
            exit_status=1, tool_log="\n".join(tool_log)))
        sys.exit(1)

    def run_tool(cmd, desc):
        log(desc)
        res = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True)
        print(res.stdout, end="")
        tool_log.append("$ " + " ".join(str(c) for c in cmd) + "\n" + res.stdout)
        return res

    def run_redirect(cmd, out_path, desc):
        log(desc)
        with open(out_path, "w") as fo:
            res = subprocess.run([str(c) for c in cmd], stdout=fo,
                                 stderr=subprocess.PIPE, text=True)
        print(res.stderr, end="")
        tool_log.append("$ " + " ".join(str(c) for c in cmd) + f" > {out_path}\n" + res.stderr)
        return res

    def out(name):
        return str(results_dir / name)

    ###########################################################################
    # Step 6: Detect compression and prepare inputs
    ###########################################################################

    log("Checking input compression...")
    r1_00 = out(f"{sample_name}_R1-00.fastq")
    kind = decompress_or_link(r1, r1_00)
    log(f"R1 input: {kind}")
    r1 = r1_00
    intermediates.append(r1)
    if not single_end:
        r2_00 = out(f"{sample_name}_R2-00.fastq")
        decompress_or_link(r2, r2_00)
        r2 = r2_00
        intermediates.append(r2)

    ###########################################################################
    # Step 7: Reformat and/or repair FASTQ (optional)
    ###########################################################################

    if reformat or repair:
        if single_end:
            r1_rf = out(f"{sample_name}_R1_reformatted-00.fastq")
            res = run_tool(["reformat.sh", f"in={r1}", f"out={r1_rf}", "tossbrokenreads=t"],
                           "Reformatting single-end FASTQ with reformat.sh...")
            if res.returncode != 0:
                fail("reformat.sh failed for R1")
            r1 = r1_rf
            intermediates.append(r1)
        else:
            r1_rf = out(f"{sample_name}_R1_reformatted-00.fastq")
            r2_rf = out(f"{sample_name}_R2_reformatted-00.fastq")
            res = run_tool(["reformat.sh", f"in={r1}", f"in2={r2}",
                            f"out={r1_rf}", f"out2={r2_rf}", "addslash=t", "tossbrokenreads=t"],
                           "Reformatting paired-end FASTQ with reformat.sh...")
            if res.returncode != 0:
                fail("reformat.sh failed for R1/R2")
            r1, r2 = r1_rf, r2_rf
            intermediates += [r1, r2]

    if repair and not single_end:
        r1_rp = out(f"{sample_name}_R1_repaired-00.fastq")
        r2_rp = out(f"{sample_name}_R2_repaired-00.fastq")
        singletons = out(f"{sample_name}_singletons_repaired-00.fastq")
        res = run_tool(["repair.sh", f"in={r1}", f"in2={r2}",
                        f"out={r1_rp}", f"out2={r2_rp}", f"outs={singletons}"],
                       "Repairing paired-end FASTQ with repair.sh...")
        if res.returncode != 0:
            fail("repair.sh failed")
        r1, r2 = r1_rp, r2_rp
        intermediates += [r1, r2, singletons]

    ###########################################################################
    # Step 8: Subsample (optional)
    ###########################################################################

    if subsample:
        r1_redu = out(f"{sample_name}_R1_redu-00.fastq")
        res = run_redirect(["seqtk", "sample", f"-s{seed}", r1, "10000"], r1_redu,
                           f"Subsampling to 10,000 reads (seed {seed})...")
        if res.returncode != 0:
            fail("seqtk subsampling R1 failed")
        r1 = r1_redu
        intermediates.append(r1)
        if not single_end:
            r2_redu = out(f"{sample_name}_R2_redu-00.fastq")
            res = run_redirect(["seqtk", "sample", f"-s{seed}", r2, "10000"], r2_redu,
                               "Subsampling R2...")
            if res.returncode != 0:
                fail("seqtk subsampling R2 failed")
            r2 = r2_redu
            intermediates.append(r2)

    ###########################################################################
    # Step 9: Adapter trimming (optional)
    ###########################################################################

    if trim_adapters:
        r1_at = out(f"{sample_name}_R1_at-01.fastq")
        if single_end:
            res = run_tool(["bbduk.sh", f"in={r1}", f"out={r1_at}", f"threads={nslots}",
                            "ktrim=r", "k=23", "mink=11", "hdist=1", f"ref={adapters_found}"],
                            "Trimming adapters (single-end)...")
        else:
            r2_at = out(f"{sample_name}_R2_at-01.fastq")
            res = run_tool(["bbduk.sh", f"in={r1}", f"in2={r2}", f"out={r1_at}", f"out2={r2_at}",
                            f"threads={nslots}", "ktrim=r", "k=23", "mink=11", "hdist=1",
                            "tpe", "tbo", f"ref={adapters_found}"],
                           "Trimming adapters (paired-end)...")
            r2 = r2_at
            intermediates.append(r2)
        if res.returncode != 0:
            fail("bbduk adapter trimming failed")
        r1 = r1_at
        intermediates.append(r1)

    ###########################################################################
    # Step 10: Quality trim paired-end reads (optional, PE only)
    ###########################################################################

    r1_qc = r2_qc = None
    if not single_end and output_pe:
        r1_qc = out(f"{sample_name}_R1_qc-02.fastq")
        r2_qc = out(f"{sample_name}_R2_qc-02.fastq")
        res = run_tool(["bbduk.sh", f"in={r1}", f"in2={r2}", f"out={r1_qc}", f"out2={r2_qc}",
                        f"minlength={min_length}", f"threads={nslots}", "qtrim=rl", f"trimq={min_qual}"],
                       "Quality trimming paired-end reads...")
        if res.returncode != 0:
            fail("bbduk quality trimming paired-end reads failed")
        if compress:
            res = run_tool(["pigz", "--keep", "--processes", nslots, r1_qc, r2_qc],
                           "Compressing paired-end QC reads...")
            if res.returncode != 0:
                fail("pigz compressing PE QC reads failed")

    ###########################################################################
    # Step 11: Merge paired-end reads (optional, PE only)
    ###########################################################################

    r_assem = r1_unassem = r2_unassem = r_discard = None
    if not single_end and output_merged:
        if merger == "pear":
            res = run_tool(["pear", "-f", r1, "-r", r2, "-o", out(sample_name),
                            "-j", nslots, "-p", pvalue, "--min-overlap", min_overlap],
                           "Merging reads with PEAR...")
            if res.returncode != 0:
                fail("Merge with PEAR failed")
            r_assem     = out(f"{sample_name}_assembled-02.fastq")
            r1_unassem  = out(f"{sample_name}_unassembled_R1-02.fastq")
            r2_unassem  = out(f"{sample_name}_unassembled_R2-02.fastq")
            r_discard   = out(f"{sample_name}_discarded-02.fastq")
            os.rename(out(f"{sample_name}.assembled.fastq"), r_assem)
            os.rename(out(f"{sample_name}.unassembled.forward.fastq"), r1_unassem)
            os.rename(out(f"{sample_name}.unassembled.reverse.fastq"), r2_unassem)
            os.rename(out(f"{sample_name}.discarded.fastq"), r_discard)
            if not os.path.getsize(r_assem):
                fail("Failed merging: 0 merged reads with PEAR")
        else:  # bbmerge
            r_assem     = out(f"{sample_name}_assembled-02.fastq")
            r1_unassem  = out(f"{sample_name}_unassembled_R1-02.fastq")
            r2_unassem  = out(f"{sample_name}_unassembled_R2-02.fastq")
            r_discard   = out(f"{sample_name}_discarded-02.fastq")
            insert_hist = out(f"{sample_name}_insert_hist.txt")
            res = run_tool(["bbmerge.sh", f"in={r1}", f"in2={r2}", f"out={r_assem}",
                            f"outu={r1_unassem}", f"outu2={r2_unassem}", f"ihist={insert_hist}"],
                           "Merging reads with BBMerge...")
            if res.returncode != 0:
                fail("Merge with BBMerge failed")
            if not os.path.getsize(r_assem):
                fail("Failed merging: 0 merged reads with BBMerge")
            Path(r_discard).touch()
        intermediates += [r_assem, r1_unassem, r2_unassem, r_discard]

    ###########################################################################
    # Step 12: Quality trim merged and unmerged reads (PE only)
    ###########################################################################

    # these files are considered intermediates, given that the final outputs are fasta files
    r_assem_qc = r1_unassem_qc = r2_unassem_qc = None
    if not single_end and output_merged:
        r_assem_qc = out(f"{sample_name}_assembled_qc-03.fastq")
        res = run_tool(["bbduk.sh", f"in={r_assem}", f"out={r_assem_qc}",
                        f"minlength={min_length}", f"threads={nslots}", "qtrim=rl", f"trimq={min_qual}"],
                       "Quality trimming merged reads...")
        if res.returncode != 0:
            fail("bbduk quality trimming merged reads failed")
        intermediates.append(r_assem_qc) 

    if not single_end and output_merged and r1_unassem and os.path.getsize(r1_unassem):
        r1_unassem_qc = out(f"{sample_name}_unassembled_R1_qc-03.fastq")
        r2_unassem_qc = out(f"{sample_name}_unassembled_R2_qc-03.fastq")
        res = run_tool(["bbduk.sh", f"in={r1_unassem}", f"in2={r2_unassem}",
                        f"out={r1_unassem_qc}", f"out2={r2_unassem_qc}",
                        f"minlength={min_length}", f"threads={nslots}", "qtrim=rl", f"trimq={min_qual}"],
                       "Quality trimming unmerged reads...")
        if res.returncode != 0:
            fail("bbduk quality trimming unmerged reads failed")
        intermediates += [r1_unassem_qc, r2_unassem_qc]

    ###########################################################################
    # Step 13: Quality trim single-end reads (SE only)
    ###########################################################################

    # fastq files as single-end reads are not removed, atlhoug there is a fasta 
    # output, because the two files are outputs of the pipeline and might be needed 
    # for downstream analysis  
    se_qc = None
    if single_end:
        se_qc = out(f"{sample_name}_se_qc-02.fastq")
        res = run_tool(["bbduk.sh", f"in={r1}", f"out={se_qc}",
                        f"minlength={min_length}", f"threads={nslots}", "qtrim=rl", f"trimq={min_qual}"],
                       "Quality trimming single-end reads...")
        if res.returncode != 0:
            fail("bbduk quality trimming single-end reads failed")
        if compress:
            res = run_tool(["pigz", "--keep", "--processes", nslots, se_qc],
                           "Compressing single-end QC reads...")
            if res.returncode != 0:
                fail("pigz compressing SE QC reads failed")

    ###########################################################################
    # Step 14: Convert final reads to FASTA
    ###########################################################################

    log("Converting to FASTA format...")

    # convert to fasta single end reads
    if single_end and se_qc and os.path.getsize(se_qc):
        se_fa = out(f"{sample_name}_se_qc-02.fasta")
        res = run_redirect(["seqtk", "seq", "-A", se_qc], se_fa, "Converting single-end reads to FASTA...")
        if res.returncode != 0:
            fail("seqtk fq2fa single-end reads failed")
        if compress:
            res = run_tool(["pigz", "--keep", "--processes", nslots, se_fa],
                           "Compressing single-end FASTA...")
            if res.returncode != 0:
                fail("pigz compressing SE FASTA failed")
            
    # convert to fasta paired-end reads
    r_assem_qc_fa = None
    if not single_end and output_merged and r_assem_qc and os.path.getsize(r_assem_qc):
        r_assem_qc_fa = out(f"{sample_name}_assembled_qc-03.fasta")
        res = run_redirect(["seqtk", "seq", "-A", r_assem_qc], r_assem_qc_fa, "Converting merged reads to FASTA...")
        if res.returncode != 0:
            fail("seqtk fq2fa merged reads failed")
        if compress:
            res = run_tool(["pigz", "--keep", "--processes", nslots, r_assem_qc_fa],
                           "Compressing merged FASTA...")
            if res.returncode != 0:
                fail("pigz compressing merged FASTA failed")

    # convert to fasta unmerged paired-end reads
    if not single_end and output_merged and r1_unassem_qc and os.path.getsize(r1_unassem_qc):
        r1_unassem_qc_fa = out(f"{sample_name}_unassembled_R1_qc-03.fasta")
        r2_unassem_qc_fa = out(f"{sample_name}_unassembled_R2_qc-03.fasta")
        res = run_redirect(["seqtk", "seq", "-A", r1_unassem_qc], r1_unassem_qc_fa, "Converting unmerged R1 to FASTA...")
        if res.returncode != 0:
            fail("seqtk fq2fa unmerged R1 failed")
        res = run_redirect(["seqtk", "seq", "-A", r2_unassem_qc], r2_unassem_qc_fa, "Converting unmerged R2 to FASTA...")
        if res.returncode != 0:
            fail("seqtk fq2fa unmerged R2 failed")
        if compress:
            res = run_tool(["pigz", "--keep", "--processes", nslots, r1_unassem_qc_fa, r2_unassem_qc_fa],
                           "Compressing unmerged FASTA...")
            if res.returncode != 0:
                fail("pigz compressing unmerged FASTA failed")

    ###########################################################################
    # Step 15: Compute stats over every intermediate FASTQ
    ###########################################################################

    log("Computing statistics...")
    rows = []
    for f in sorted(glob.glob(str(results_dir / "*.fastq")) + glob.glob(str(results_dir / "*.fastq.gz"))):
        base = os.path.basename(f)
        n = count_fastq(f)
        length = mean_length(f, fmt="fastq") if os.path.getsize(f) else 0
        order = parse_order(base)
        rows.append({"sample": sample_name, "file": base, "stat": "num_seq",     "value": n,      "order": order})
        rows.append({"sample": sample_name, "file": base, "stat": "mean_length", "value": length, "order": order})

    with open(stats_out, "w") as fh:
        fh.write("sample\tfile\tstat\tvalue\torder\n")
        for r in rows:
            fh.write(f"{r['sample']}\t{r['file']}\t{r['stat']}\t{r['value']}\t{r['order']}\n")

    ###########################################################################
    # Step 16: Remove the plain QC-trimmed reads once compressed
    ###########################################################################

    # pigz --keep leaves the plain file in place alongside the new .gz; keeping
    # both around is redundant once compressed. 
    # Deferred to here (rather than done inline right after each pigz
    # call above) because Step 14's FASTA conversion still needs the plain se_qc; 
    # r1_qc/r2_qc (when using output merged) aren't read again after Step 14, 
    # so removing them here too is just for a single, easy-to-audit cleanup point.

    if compress:
        # remove uncompressed fastq and fasta single-end reads 
        if single_end:
            if se_qc and os.path.exists(se_qc):
                os.remove(se_qc)
            if se_fa and os.path.exists(se_fa):
                os.remove(se_fa)

        # remove uncompressed fastq and fasta merged paired-reads
        if not single_end and output_merged:
            if r_assem_qc and os.path.exists(r_assem_qc):
                os.remove(r_assem_qc)
            if r_assem_qc_fa and os.path.exists(r_assem_qc_fa):
                os.remove(r_assem_qc_fa)

        # remove uncompressed fastq and fasta unmerged reads
        if not single_end and output_merged:
            if r1_unassem_qc and os.path.exists(r1_unassem_qc):
                os.remove(r1_unassem_qc)
            if r1_unassem_qc_fa and os.path.exists(r1_unassem_qc_fa):
                os.remove(r1_unassem_qc_fa)
            if r2_unassem_qc and os.path.exists(r2_unassem_qc):
                os.remove(r2_unassem_qc)
            if r2_unassem_qc_fa and os.path.exists(r2_unassem_qc_fa):
                os.remove(r2_unassem_qc_fa)

        # remove uncompressed paired reads
        if not single_end and output_pe:
            if r1_qc and os.path.exists(r1_qc):
                os.remove(r1_qc)
            if r2_qc and os.path.exists(r2_qc):
                os.remove(r2_qc)

    ###########################################################################
    # Step 17: Clean intermediates (optional)
    ###########################################################################

    # Files that must survive cleanup: the emitted QC reads for downstream
    # assembly (PE: *_qc-02.fastq[.gz]; SE: *_se_qc-02.fastq[.gz]) and the FASTA
    # that becomes <sample>_workable.
    if clean:
        log("Cleaning intermediate files...")
        for f in intermediates:
            if f and os.path.realpath(f):
                os.remove(f)

    ###########################################################################
    # Step 18: Write log
    ###########################################################################

    gz_suffix = ".gz" if compress else ""
    outputs = [f"Statistics: {stats_out}"]
    if  output_merged and r_assem_qc_fa:
        outputs.append(f"Merged FASTA: {r_assem_qc_fa}{gz_suffix}")
    if  output_merged and r1_unassem_qc_fa:
        outputs.append(f"Unmerged FASTA R1: {r1_unassem_qc_fa}{gz_suffix}")
        outputs.append(f"Unmerged FASTA R2: {r2_unassem_qc_fa}{gz_suffix}")
    if not single_end and output_pe and r1_qc:
        outputs.append(f"QC R1 FASTQ: {r1_qc}{gz_suffix}")
        outputs.append(f"QC R2 FASTQ: {r2_qc}{gz_suffix}")
    if single_end and se_qc:
        outputs.append(f"QC single-end reads FASTQ: {se_qc}{gz_suffix}")
        outputs.append(f"QC single-end reads FASTA: {se_fa}{gz_suffix}")

    log("\033[0;32m2-preprocess.py completed successfully\033[0m")
    log_out.write_text(build_log(
        script_name=SCRIPT_NAME,
        script_desc=SCRIPT_DESC,
        sample_name=sample_name,
        inputs=inputs,
        params=params,
        outputs=outputs, command=command, exit_status=0,
        tool_log="\n".join(tool_log)))

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
