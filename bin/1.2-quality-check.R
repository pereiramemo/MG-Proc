#!/usr/bin/env Rscript

###############################################################################
### 1. Set env
###############################################################################

suppressMessages({
  library(tidyverse)
  library(ShortRead)
  library(doParallel)
  library(dada2)
  library(optparse)
  library(this.path)
})

# Import shared helpers from bin/utils.R (sibling module).
source(file.path(dirname(this.path::this.path()), "utils.R"))

script_name <- "1.2-quality-check.R"
script_desc <- paste(
  "Comparative QC plots across all samples (mean quality vs read count,",
  "count histograms, PhiX contamination)."
)

###############################################################################
### 2. Parse command line arguments
###############################################################################

option_list <- list(
  make_option("--input_dir",
    type = "character", default = NULL,
    help = paste(
      "Directory containing the FASTQ files listed in --input_tsv",
      "(files are located by basename)"
    ),
    metavar = "character"
  ),
  make_option("--input_tsv",
    type = "character", default = NULL,
    help = paste(
      "TSV samplesheet with columns sample_name, reads1, reads2",
      "(tab-delimited). Only the basename of reads1/reads2 is used to locate",
      "files under --input_dir. Leave reads2 empty for single-end samples;",
      "a sheet must be either all paired-end or all single-end."
    ),
    metavar = "character"
  ),
  make_option("--output_dir",
    type = "character", default = NULL,
    help = "Output directory", metavar = "character"
  ),
  make_option("--nslots",
    type = "integer", default = 12,
    help = "Number of threads to use [default=%default]",
    metavar = "integer"
  ),
  make_option("--sample_size",
    type = "integer", default = 10000,
    help = paste(
      "Reads to subsample per file for quality and PhiX",
      "estimation [default=%default]"
    ),
    metavar = "integer"
  ),
  make_option("--overwrite",
    type = "character", default = "f",
    help = "Overwrite previous output [default=%default]",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_bool <- function(x, flag) {
  if (tolower(x) %in% c("t", "true")) {
    return(TRUE)
  }
  if (tolower(x) %in% c("f", "false")) {
    return(FALSE)
  }
  stop(paste0("--", flag, " must be t/f (got '", x, "')"), call. = FALSE)
}
overwrite <- parse_bool(opt$overwrite, "overwrite")

if (is.null(opt$input_dir) || is.null(opt$input_tsv) || is.null(opt$output_dir)) {
  log_error("--input_dir, --input_tsv, and --output_dir are required arguments.")
  quit(status = 1)
}

input_dir <- opt$input_dir
orig_input_dir <- input_dir # kept for the log's "Input data" section; input_dir
                            # itself is repointed to a decompressed staging dir below
input_tsv <- opt$input_tsv
output_dir <- opt$output_dir
nslots <- opt$nslots
sample_size <- opt$sample_size

# Dev only section
# input_dir <- "/home/epereira/workspace/repos/tools/MG-Proc/tests/data_samo" # nolintr
# output_dir <- "/home/epereira/workspace/repos/tools/MG-Proc/tests/output/1.2-quality-check-out" # nolintr
# input_tsv <- "/home/epereira/workspace/repos/tools/MG-Proc/tests/data_samo/samplesheet.tsv" # nolintr
# nslots <- 12
# sample_size <- 10000

###############################################################################
### 3. Validate input directory, samplesheet, and dependencies
###############################################################################

if (!dir.exists(input_dir)) {
  log_error(paste("Input directory does not exist:", input_dir))
  quit(status = 1)
}

if (!file.exists(input_tsv)) {
  log_error(paste("Input TSV does not exist:", input_tsv))
  quit(status = 1)
}

samplesheet <- read_tsv(input_tsv, col_types = cols(.default = "c"), progress = FALSE)
required_cols <- c("sample_name", "reads1", "reads2")
if (!all(required_cols %in% names(samplesheet))) {
  log_error(paste0(
    "--input_tsv must have columns: ", paste(required_cols, collapse = ", "),
    " (got: ", paste(names(samplesheet), collapse = ", "), ")"
  ))
  quit(status = 1)
}
if (nrow(samplesheet) == 0) {
  log_error(paste("--input_tsv has no data rows:", input_tsv))
  quit(status = 1)
}
if (anyDuplicated(samplesheet$sample_name)) {
  log_error("--input_tsv has duplicate sample_name values.")
  quit(status = 1)
}

samplesheet$reads2[is.na(samplesheet$reads2)] <- ""
single_end_flags_list <- trimws(samplesheet$reads2) == ""
if (length(unique(single_end_flags_list)) > 1) {
  log_error(paste(
    "--input_tsv mixes single-end and paired-end rows",
    "(reads2 must be either always empty or always populated)."
  ))
  quit(status = 1)
}
single_end_flag <- single_end_flags_list[1]

registerDoParallel(cores = nslots)

###############################################################################
### 4. Prepare output directories (output/, logs/, stats/)
###############################################################################

if (dir.exists(output_dir)) {
  if (!overwrite) {
    log_error(paste(
      "Output directory already exists:", output_dir,
      "- use --overwrite t to overwrite"
    ))
    quit(status = 1)
  }
  log_warn(paste("Overwriting existing directory:", output_dir))
  unlink(output_dir, recursive = TRUE)
}

results_dir <- file.path(output_dir, "output")
logs_dir <- file.path(output_dir, "logs")
stats_dir <- file.path(output_dir, "stats")
for (d in c(results_dir, logs_dir, stats_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_out <- file.path(logs_dir, "1.2-quality-check.log")
stats_out <- file.path(stats_dir, "1.2-quality-check-stats.tsv")

###############################################################################
### 5. Find and stage input files
###############################################################################

locate_read_file <- function(path, sample_name, mate_label) {
  f <- file.path(input_dir, basename(path))
  if (!file.exists(f)) {
    log_error(paste0(
      "Sample '", sample_name, "': ", mate_label, " file not found: ",
      basename(path), " (looked in ", input_dir, ")"
    ))
    quit(status = 1)
  }
  f
}

samplesheet$reads1_path <- mapply(
  locate_read_file, samplesheet$reads1, samplesheet$sample_name,
  MoreArgs = list(mate_label = "reads1")
)
if (!single_end_flag) {
  samplesheet$reads2_path <- mapply(
    locate_read_file, samplesheet$reads2, samplesheet$sample_name,
    MoreArgs = list(mate_label = "reads2")
  )
}

log_msg(if (single_end) {
  paste("Found", nrow(samplesheet), "single-end files")
} else {
  paste("Found", nrow(samplesheet), "samples (R1 + R2 each)")
})

# ShortRead's qa()/countFastq()/FastqSampler() only understand gzip (or
# plain), not bzip2, so stage a plain/decompressed copy of every file, named
# after its sample (<sample_name>_R1.fastq[/_R2.fastq]) rather than its
# original filename. decompress_or_link() picks gzip/bzip2/plain by sniffing
# the file's magic bytes, so any destination name works. Because the staged
# names are our own convention rather than something derived from user input,
# pattern_r1/pattern_r2 below are fixed constants - only qa()'s dirPath+pattern
# API (see ShortRead docs) still needs a "pattern" at all.
pattern_r1 <- "_R1\\.fastq$"
pattern_r2 <- "_R2\\.fastq$"

staged_dir <- tempfile(pattern = "1.2-quality-check-", tmpdir = output_dir)
dir.create(staged_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(samplesheet))) {
  decompress_or_link(
    samplesheet$reads1_path[i],
    file.path(staged_dir, paste0(samplesheet$sample_name[i], "_R1.fastq"))
  )
  if (!single_end) {
    decompress_or_link(
      samplesheet$reads2_path[i],
      file.path(staged_dir, paste0(samplesheet$sample_name[i], "_R2.fastq"))
    )
  }
}

input_dir <- staged_dir

raw_r1 <- sort(file.path(staged_dir, paste0(samplesheet$sample_name, "_R1.fastq")))
if (!single_end) {
  raw_r2 <- sort(file.path(staged_dir, paste0(samplesheet$sample_name, "_R2.fastq")))
}

###############################################################################
### 6. Per-sample read counts (R1 / single-end)
###############################################################################

seq_counts_df <- foreach(i = raw_r1, .combine = rbind) %dopar% {
  count_seqs(i, pattern_r1)
}

###############################################################################
### 7. R1 mean quality vs read count
###############################################################################

x_r1 <- qa(
  dirPath = input_dir, pattern = pattern_r1,
  sample = TRUE, n = sample_size
)
qa_means <- x_r1[["perCycle"]][["quality"]] |>
  group_by(lane) |>
  summarize(mean_q = sum(Score * Count) / sum(Count), .groups = "drop")
qa_means$lane <- sapply(
  qa_means$lane,
  function(x) extract_sample_name(x, pattern_r1)
)
qa_means2counts <- left_join(qa_means, seq_counts_df,
                             by = c("lane" = "sample"))

text_size <- 2
p_r1 <- ggplot(qa_means2counts, aes(x = mean_q, y = nseq)) +
  geom_point() +
  scale_y_log10() +
  ylab("Read counts (log)") +
  xlab("Mean quality score (R1)") +
  geom_text(aes(label = as.character(lane)),
    hjust = 0.5, vjust = -1, size = text_size
  )

ggsave(p_r1,
  filename = file.path(results_dir, "r1_mean_q_vs_nseq.png"),
  device = "png", width = 5, height = 4, dpi = 300
)

###############################################################################
### 8. R2 mean quality vs read count (paired-end only)
###############################################################################

if (!single_end) {
  x_r2 <- qa(
    dirPath = input_dir, pattern = pattern_r2,
    sample = TRUE, n = sample_size
  )
  qa_means_r2 <- x_r2[["perCycle"]][["quality"]] |>
    group_by(lane) |>
    summarize(mean_q = sum(Score * Count) / sum(Count), .groups = "drop")
  qa_means_r2$lane <- sapply(
    qa_means_r2$lane,
    function(x) extract_sample_name(x, pattern_r2)
  )
  qa_means2counts_r2 <- left_join(qa_means_r2, seq_counts_df,
                                  by = c("lane" = "sample"))

  p_r2 <- ggplot(qa_means2counts_r2, aes(x = mean_q, y = nseq)) +
    geom_point() +
    scale_y_log10() +
    ylab("Read counts (log)") +
    xlab("Mean quality score (R2)") +
    geom_text(aes(label = as.character(lane)),
      hjust = 0.5, vjust = -1, size = text_size
    )

  ggsave(p_r2,
    filename = file.path(results_dir, "r2_mean_q_vs_nseq.png"),
    device = "png", width = 5, height = 4, dpi = 300
  )
}

###############################################################################
### 9. Read-count histograms
###############################################################################

samples_hist_p <- ggplot(qa_means2counts, aes(nseq)) +
  geom_histogram(bins = 40) +
  ylab("Num of samples")
samples_hist_log_p <- ggplot(qa_means2counts, aes(log(nseq))) +
  geom_histogram(bins = 40) +
  ylab("Num of samples")

ggsave(samples_hist_p,
  filename = file.path(results_dir, "samples_hist.png"),
  device = "png", width = 5, height = 4, dpi = 300
)
ggsave(samples_hist_log_p,
  filename = file.path(results_dir, "samples_hist_log.png"),
  device = "png", width = 5, height = 4, dpi = 300
)

###############################################################################
### 10. PhiX contamination estimate
###############################################################################

phix_counts_df <- foreach(i = raw_r1, .combine = rbind) %dopar% {
  count_phix_seqs(i, n_sample = sample_size, seed = 123)
}

phix_df <- data.frame(
  sample = sapply(
    phix_counts_df$file,
    function(x) extract_sample_name(x, pattern_r1)
  ),
  perc = 100 * phix_counts_df$nphix / phix_counts_df$n,
  nphix = phix_counts_df$nphix
)

phix_plot <- phix_df
undet <- grep(pattern = "Undetermined", phix_plot$sample)
if (length(undet) > 0) {
  phix_plot <- phix_plot[-undet, ]
}
phix_plot$color <- ifelse(phix_plot$perc > 0.005, "indianred", "gray50")

perc_phix_barplot <- ggplot(
  phix_plot,
  aes(x = sample, y = perc, fill = color)
) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(
    values = c("gray50" = "gray50", "indianred" = "indianred"),
    labels = c("gray50" = "Phix <= 0.005%", "indianred" = "Phix > 0.005%"),
    name = ""
  ) +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(perc_phix_barplot,
  filename = file.path(results_dir, "samples_perc_phix_barplot.png"),
  device = "png", width = 10, height = 4, dpi = 300
)

###############################################################################
### 11. Stats table (samples as rows, statistics as columns)
###############################################################################

stats_tbl <- qa_means2counts |>
  dplyr::rename(sample = lane, mean_q_r1 = mean_q)
if (!single_end) {
  stats_tbl <- stats_tbl |>
    left_join(
      qa_means2counts_r2 |>
        dplyr::rename(sample = lane, mean_q_r2 = mean_q) |>
        select(sample, mean_q_r2),
      by = "sample"
    )
}
phix_summary <- phix_df |> transmute(sample = sample, phix_pct = perc)
stats_tbl <- stats_tbl |>
  left_join(phix_summary, by = "sample") |>
  dplyr::select(sample, nseq, mean_q_r1, dplyr::any_of("mean_q_r2"), phix_pct)

write_tsv(stats_tbl, stats_out)

unlink(staged_dir, recursive = TRUE)

###############################################################################
### 12. Log and summary
###############################################################################

log_msg(paste("Processed", length(raw_r1), "samples"))
log_msg("\033[0;32m1.2-quality-check.R completed successfully\033[0m")

generated <- c(
  "r1_mean_q_vs_nseq.png",
  if (!single_end) "r2_mean_q_vs_nseq.png",
  "samples_hist.png", "samples_hist_log.png",
  "samples_perc_phix_barplot.png"
)

inputs <- c(
  paste("Input directory:", orig_input_dir),
  paste("Input samplesheet:", input_tsv),
  paste("Samples:", paste(samplesheet$sample_name, collapse = ", "))
)

params <- c(
  paste("Threads:", nslots),
  paste("Read type:", if (single_end) "single-end" else "paired-end"),
  paste("Sample size:", sample_size)
)

outputs <- c(
  file.path(results_dir, generated),
  paste("Statistics:", stats_out)
)
command <- paste(
  c(script_name, commandArgs(trailingOnly = TRUE)),
  collapse = " "
)

log_lines <- build_log(
  script_name, script_desc,
  sample_name = "all samples",
  inputs = inputs, params = params, outputs = outputs,
  command = command, exit_status = 0,
  tool_log = paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(log_lines, log_out)
