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
    help = "Input directory with FASTQ files", metavar = "character"
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
  make_option("--single_end",
    type = "character", default = "f",
    help = paste(
      "Process single-end reads instead of paired-end",
      "[default=%default]"
    ),
    metavar = "character"
  ),
  make_option("--reads_pattern",
    type = "character", default = "*_{1,2}.fastq.gz",
    help = paste(
      "Paired-end glob; use a {1,2} token to mark the R1/R2 mates",
      "(e.g. *_{1,2}.fastq.gz). Ignored when --single_end t [default=%default]"
    ),
    metavar = "character"
  ),
  make_option("--se_reads_pattern",
    type = "character", default = "*.fastq.gz",
    help = paste(
      "Single-end glob matching the read files directly (e.g. *.fastq.gz).",
      "Used only when --single_end t [default=%default]"
    ),
    metavar = "character"
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
single_end <- parse_bool(opt$single_end, "single_end")
overwrite <- parse_bool(opt$overwrite, "overwrite")

if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  log_error("--input_dir and --output_dir are required arguments.")
  quit(status = 1)
}

input_dir <- opt$input_dir
orig_input_dir <- input_dir # kept for the log's "Input data" section; input_dir
                            # itself is repointed to a decompressed staging dir below
output_dir <- opt$output_dir
nslots <- opt$nslots
sample_size <- opt$sample_size
reads_pattern <- opt$reads_pattern
se_reads_pattern <- opt$se_reads_pattern

# Derive the read-identifying suffix(es) from the input glob(s). The leading "*"
# is the sample name; the remainder identifies the read. In paired-end mode a
# {1,2} token in --reads_pattern marks the R1/R2 mate position; in single-end
# mode --se_reads_pattern matches the read files directly. The suffixes are used
# as substring patterns by list.files()/qa() and to strip the sample name.
if (single_end) {
  pattern_r1 <- sub("^\\*", "", se_reads_pattern)
  pattern_r2 <- NULL
} else {
  read_suffix <- sub("^\\*", "", reads_pattern)
  if (!grepl("{1,2}", read_suffix, fixed = TRUE)) {
    log_error(paste0(
      "--reads_pattern must contain a {1,2} token in paired-end mode (got '",
      reads_pattern, "')"
    ))
    quit(status = 1)
  }
  pattern_r1 <- sub("{1,2}", "1", read_suffix, fixed = TRUE)
  pattern_r2 <- sub("{1,2}", "2", read_suffix, fixed = TRUE)
}

# Dev only section
# input_dir <- "/home/epereira/workspace/repos/tools/MG-Proc/tests/data/SRR12479690/" # nolintr
# output_dir <- "/home/epereira/workspace/repos/tools/MG-Proc/tests/output/1.2-quality-check-out" # nolintr
# nslots <- 12
# single_end <- FALSE
# pattern_r1 <- "_1.fastq.gz"
# pattern_r2 <- "_2.fastq.gz"
# sample_size <- 10000

###############################################################################
### 3. Validate input directory and dependencies
###############################################################################

if (!dir.exists(input_dir)) {
  log_error(paste("Input directory does not exist:", input_dir))
  quit(status = 1)
}

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

raw_r1 <- sort(list.files(input_dir, pattern = pattern_r1, full.names = TRUE))
if (length(raw_r1) == 0) {
  log_error(paste("No files found matching pattern:", pattern_r1))
  quit(status = 1)
}

if (!single_end) {
  raw_r2 <- sort(
    list.files(input_dir, pattern = pattern_r2, full.names = TRUE)
  )
  if (length(raw_r2) == 0) {
    log_error(paste("No R2 files found matching pattern:", pattern_r2))
    quit(status = 1)
  }
  if (length(raw_r1) != length(raw_r2)) {
    log_error(paste(
      "Mismatch in number of files. R1:", length(raw_r1),
      "R2:", length(raw_r2)
    ))
    quit(status = 1)
  }
  log_msg(paste(
    "Found", length(raw_r1), "R1 files and", length(raw_r2), "R2 files"
  ))
} else {
  log_msg(paste("Found", length(raw_r1), "single-end files"))
}

# ShortRead's qa()/countFastq()/FastqSampler() only understand gzip (or
# plain), not bzip2, so stage a plain/decompressed copy of every matched file
# and point the rest of the script at that staging dir instead of input_dir.
# Compression suffixes are stripped from the (derived) patterns to match the
# staged, decompressed filenames.
strip_compression_ext <- function(x) sub("\\.(gz|bz2)$", "", x)

staged_dir <- tempfile(pattern = "1.2-quality-check-", tmpdir = output_dir)
dir.create(staged_dir, recursive = TRUE, showWarnings = FALSE)

for (f in raw_r1) {
  decompress_or_link(f, file.path(staged_dir, strip_compression_ext(basename(f))))
}
if (!single_end) {
  for (f in raw_r2) {
    decompress_or_link(f, file.path(staged_dir, strip_compression_ext(basename(f))))
  }
}

orig_pattern_r1 <- pattern_r1
orig_pattern_r2 <- pattern_r2
pattern_r1 <- strip_compression_ext(pattern_r1)
if (!single_end) pattern_r2 <- strip_compression_ext(pattern_r2)
input_dir <- staged_dir

raw_r1 <- sort(list.files(input_dir, pattern = pattern_r1, full.names = TRUE))
if (!single_end) {
  raw_r2 <- sort(list.files(input_dir, pattern = pattern_r2, full.names = TRUE))
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
  paste("Reads pattern:", if (single_end) se_reads_pattern else reads_pattern),
  paste("R1 pattern (derived):", orig_pattern_r1),
  if (!single_end) paste("R2 pattern (derived):", orig_pattern_r2)
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
