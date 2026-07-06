#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import glob
import gzip
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# Import shared helpers from bin/toolbox.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolbox import log, log_warn, log_error, derive_sample_name, build_log, count_fasta

SCRIPT_NAME = "3-assembly-and-map.py"
SCRIPT_DESC = ("De novo assembly (MEGAHIT) and read mapping (BWA-MEM + SAMtools) for a "
               "single metagenomic sample, with optional Picard duplicate removal.")

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(description=SCRIPT_DESC)
    p.add_argument("--reads1",            required=True, help="R1 (or single-end) reads, fastq/fa (required)")
    p.add_argument("--reads2",            default=None,  help="R2 reads (paired-end only)")
    p.add_argument("--single_end",        choices=["t", "f"], default="f", help="Process single-end reads [default=f]")
    p.add_argument("--sample_name",       default="metagenomex", help="Sample name for output files [default=metagenomex]")
    p.add_argument("--contigs",           default=None,  help="Pre-assembled contigs FASTA (takes precedence over --assem_dir)")
    p.add_argument("--assem_dir",         default=None,  help="Directory with previously computed assemblies")
    p.add_argument("--assem_preset",      default="meta-sensitive", help="MEGAHIT preset [default=meta-sensitive]")
    p.add_argument("--nslots",            type=int, default=12,  help="Number of threads [default=12]")
    p.add_argument("--min_contig_length", type=int, default=250, help="Minimum contig length to keep [default=250]")
    p.add_argument("--output_dir",        default="mg-clust_output-1", help="Output directory [default=mg-clust_output-1]")
    p.add_argument("--overwrite",         choices=["t", "f"], default="f", help="Overwrite previous folder [default=f]")
    p.add_argument("--remove_duplicates", choices=["t", "f"], default="f", help="Remove PCR duplicates with Picard [default=f]")
    return p.parse_args()


def find_contigs_file(assem_dir, sample_name):
    """Reproduce the bash find_contigs_file search (sample-prefixed names first,
    then any *.contigs.{fa,fasta,fna}[.gz])."""
    extensions = ["contigs.fa", "contigs.fasta", "contigs.fna", "fa", "fasta", "fna",
                  "contigs.fa.gz", "contigs.fasta.gz", "contigs.fna.gz",
                  "fa.gz", "fasta.gz", "fna.gz"]
    search_paths = [os.path.join(assem_dir, sample_name), assem_dir]
    for d in search_paths:
        for ext in extensions:
            cand = os.path.join(d, f"{sample_name}.{ext}")
            if os.path.isfile(cand):
                return cand
    for d in search_paths:
        for ext in ("fa", "fasta", "fna"):
            for pattern in (f"*.contigs.{ext}", f"*.contigs.{ext}.gz"):
                hits = sorted(glob.glob(os.path.join(d, pattern)))
                if hits:
                    log_warn(f"Found contigs file '{hits[0]}' not matching '{sample_name}.contigs.*'; using it anyway.")
                    return hits[0]
    return None


################################################################################
# 3. Define main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    r1                = opts.reads1
    r2                = opts.reads2
    single_end        = opts.single_end == "t"
    sample_name       = opts.sample_name
    contigs           = opts.contigs
    assem_dir         = opts.assem_dir
    assem_preset      = opts.assem_preset
    nslots            = opts.nslots
    min_contig_length = opts.min_contig_length
    output_dir        = Path(opts.output_dir)
    overwrite         = opts.overwrite == "t"
    remove_duplicates = opts.remove_duplicates == "t"

    tool_log = []

    ###########################################################################
    # Step 2: Validate inputs, tools, and prepare output directories
    ###########################################################################

    if not os.path.isfile(r1):
        log_error(f"--reads1 file does not exist: {r1}")
        sys.exit(1)
    if not single_end:
        if not r2:
            log_error("--reads2 is required for paired-end mode.")
            sys.exit(1)
        if not os.path.isfile(r2):
            log_error(f"--reads2 file does not exist: {r2}")
            sys.exit(1)

    log("Checking dependencies...")
    required = ["bwa", "samtools"]
    if not contigs and not assem_dir:
        required.append("megahit")
    if remove_duplicates:
        required.append("picard")
    missing = [t for t in required if not shutil.which(t)]
    if missing:
        log_error("Missing required tools: " + ", ".join(missing))
        sys.exit(1)

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

    log_out   = logs_dir  / f"3-assembly-and-map-{sample_name}.log"
    stats_out = stats_dir / f"3-assembly-and-map-{sample_name}-stats.tsv"

    mode = "single-end" if single_end else "paired-end"
    inputs = [f"R1: {r1}"] + ([] if single_end else [f"R2: {r2}"])
    params = [
        f"Sample name: {sample_name}", f"Read type: {mode}",
        f"Assembly preset: {assem_preset}", f"Min contig length: {min_contig_length}",
        f"Threads: {nslots}", f"Remove duplicates: {'yes' if remove_duplicates else 'no'}",
        f"External contigs: {contigs or assem_dir or 'no (de novo assembly)'}",
    ]
    command = " ".join([SCRIPT_NAME] + sys.argv[1:])

    def write_log(exit_status, outputs):
        log_out.write_text(build_log(
            script_name=SCRIPT_NAME,
            script_desc=SCRIPT_DESC,
            sample_name=sample_name,
            inputs=inputs,
            params=params,
            outputs=outputs, command=command, exit_status=exit_status,
            tool_log="\n".join(tool_log)))

    def fail(msg):
        log_error(msg)
        write_log(1, ["(failed before completion)"])
        sys.exit(1)

    def run_tool(cmd, desc):
        log(desc)
        res = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True)
        print(res.stdout, end="")
        tool_log.append("$ " + " ".join(str(c) for c in cmd) + "\n" + res.stdout)
        return res

    def run_shell(cmd_str, desc):
        log(desc)
        res = subprocess.run(cmd_str, shell=True, executable="/bin/bash",
                             stderr=subprocess.PIPE, text=True)
        print(res.stderr, end="")
        tool_log.append("$ " + cmd_str + "\n" + res.stderr)
        return res

    ###########################################################################
    # Step 3: De novo assembly or resolve external contigs
    ###########################################################################

    assembly_file = None

    if contigs:
        log(f"Using externally provided contigs file: {contigs}")
        if not os.path.isfile(contigs):
            fail(f"Contigs file does not exist: {contigs}")
        assembly_file = contigs
    elif assem_dir:
        log(f"Searching for contigs in assembly directory: {assem_dir}")
        if not os.path.isdir(assem_dir):
            fail(f"Assembly directory does not exist: {assem_dir}")
        assembly_file = find_contigs_file(assem_dir, sample_name)
        if not assembly_file:
            fail(f"Could not find contigs file for sample '{sample_name}' in '{assem_dir}'. "
                 "Use --contigs to specify the exact path.")
        log(f"Found contigs file: {assembly_file}")
    else:
        log("Running MEGAHIT assembly...")
        megahit_dir = output_dir / "megahit_tmp"
        reads = ["-r", r1] if single_end else ["-1", r1, "-2", r2]
        cmd = ["megahit", "--num-cpu-threads", nslots, "--presets", assem_preset,
               "--min-contig-len", min_contig_length, "--out-prefix", sample_name,
               "--out-dir", str(megahit_dir)] + reads
        res = run_tool(cmd, "Assembling with MEGAHIT...")
        if res.returncode != 0:
            fail("MEGAHIT assembly failed.")
        produced = megahit_dir / f"{sample_name}.contigs.fa"
        if not os.path.isfile(produced):
            fail("MEGAHIT did not produce a contigs file.")
        assembly_file = str(results_dir / f"{sample_name}.contigs.fa")
        shutil.copyfile(produced, assembly_file)
        shutil.rmtree(megahit_dir, ignore_errors=True)

    # Copy/decompress an external assembly into output/ so BWA indices are local.
    if contigs or assem_dir:
        if assembly_file.endswith(".gz"):
            local = str(results_dir / os.path.basename(assembly_file)[:-3])
            log("Decompressing external assembly into output directory...")
            with gzip.open(assembly_file, "rb") as fi, open(local, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        else:
            local = str(results_dir / os.path.basename(assembly_file))
            shutil.copyfile(assembly_file, local)
        assembly_file = local

    ###########################################################################
    # Step 4: Guard on contig count
    ###########################################################################

    n_contigs = count_fasta(assembly_file)
    log(f"Assembled contigs: {n_contigs}")
    if n_contigs < 5:
        log_warn(f"Not enough assembled sequences to continue ({n_contigs} < 5). Exiting gracefully.")
        stats_out.write_text("sample\tn_contigs\tn_mapped_reads\n"
                             f"{sample_name}\t{n_contigs}\t0\n")
        write_log(0, [f"Contigs: {assembly_file}", f"Statistics: {stats_out}"])
        sys.exit(0)

    ###########################################################################
    # Step 5: Map reads (BWA-MEM -> BAM, sort, index)
    ###########################################################################

    res = run_tool(["bwa", "index", assembly_file], "Indexing assembly with BWA...")
    if res.returncode != 0:
        fail("bwa index failed.")

    bam       = str(results_dir / f"{sample_name}.bam")
    sorted_bam = str(results_dir / f"{sample_name}_sorted.bam")
    q_asm = shlex.quote(assembly_file)
    reads_str = shlex.quote(r1) if single_end else f"{shlex.quote(r1)} {shlex.quote(r2)}"
    res = run_shell(
        f"bwa mem -M -t {nslots} {q_asm} {reads_str} | "
        f"samtools view -@ {nslots} -q 10 -F 260 -b > {shlex.quote(bam)}",
        "Mapping reads with BWA-MEM and filtering (q>=10, primary alignments)...")
    if res.returncode != 0:
        fail("bwa mem or samtools view failed.")

    res = run_tool(["samtools", "sort", "-@", nslots, "-o", sorted_bam, bam], "Sorting BAM...")
    if res.returncode != 0:
        fail("samtools sort failed.")
    res = run_tool(["samtools", "index", "-@", nslots, sorted_bam], "Indexing sorted BAM...")
    if res.returncode != 0:
        fail("samtools index failed.")

    ###########################################################################
    # Step 6: Remove duplicates with Picard (optional)
    ###########################################################################

    if remove_duplicates:
        markdup_bam = str(results_dir / f"{sample_name}_sorted_markdup.bam")
        metrics     = str(results_dir / f"{sample_name}_sorted_markdup.metrics.txt")
        tmp_dir     = output_dir / "tmp"
        tmp_dir.mkdir(exist_ok=True)
        res = run_tool([
            "picard", "MarkDuplicates",
            f"INPUT={sorted_bam}", f"OUTPUT={markdup_bam}", f"METRICS_FILE={metrics}",
            "REMOVE_DUPLICATES=TRUE", "ASSUME_SORTED=TRUE",
            "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900", f"TMP_DIR={tmp_dir}",
        ], "Marking and removing duplicates with Picard...")
        if res.returncode != 0:
            fail("Picard MarkDuplicates failed.")
        res = run_tool(["samtools", "index", "-@", nslots, markdup_bam], "Indexing duplicate-removed BAM...")
        if res.returncode != 0:
            fail("samtools index on markdup BAM failed.")
        final_bam = markdup_bam
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        log("Skipping duplicate removal (--remove_duplicates f)")
        final_bam = sorted_bam

    ###########################################################################
    # Step 7: Stats
    ###########################################################################

    res = subprocess.run(["samtools", "view", "-c", final_bam],
                         stdout=subprocess.PIPE, text=True)
    n_mapped = res.stdout.strip() if res.returncode == 0 else "NA"
    stats_out.write_text("sample\tn_contigs\tn_mapped_reads\n"
                        f"{sample_name}\t{n_contigs}\t{n_mapped}\n")

    ###########################################################################
    # Step 8: Clean intermediates
    ###########################################################################

    log("Removing intermediate mapping files...")
    if os.path.exists(bam):
        os.remove(bam)
    if remove_duplicates:
        for f in (sorted_bam, sorted_bam + ".bai"):
            if os.path.exists(f):
                os.remove(f)
    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        idx = assembly_file + ext
        if os.path.exists(idx):
            os.remove(idx)

    ###########################################################################
    # Step 9: Write log
    ###########################################################################

    outputs = [f"Contigs: {assembly_file}", f"Final BAM: {final_bam}",
               f"Statistics: {stats_out}"]
    log(f"\033[0;32m3-assembly-and-map.py completed successfully\033[0m (final BAM: {final_bam})")
    write_log(0, outputs)


################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
