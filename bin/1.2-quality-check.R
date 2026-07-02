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

# Import shared helpers from bin/toolbox.R (sibling module).
source(file.path(dirname(this.path::this.path()), "toolbox.R"))

SCRIPT_NAME <- "1.2-quality-check.R"
SCRIPT_DESC <- "Comparative QC plots across all samples (mean quality vs read count, count histograms, PhiX contamination)." # nolintr

###############################################################################
### 2. Parse command line arguments
###############################################################################

option_list <- list(
  make_option(c("--input_dir"),  type = "character", default = NULL,
              help = "Input directory with FASTQ files", metavar = "character"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory", metavar = "character"),
  make_option(c("--nslots"),     type = "integer", default = 12,
              help = "Number of threads to use [default=%default]", metavar = "integer"), 
  make_option(c("--single_end"), type = "character", default = "f",
              help = "Process single-end reads instead of paired-end [default=%default]", metavar = "character"),
  make_option(c("--r1_pattern"), type = "character", default = "_R1_001.fastq.gz",
              help = "Pattern for R1 FASTQ files, or single-end files when --single_end t [default=%default]", metavar = "character"),
  make_option(c("--r2_pattern"), type = "character", default = "_R2_001.fastq.gz",
              help = "Pattern for R2 FASTQ files (ignored when --single_end t) [default=%default]", metavar = "character"),
  make_option(c("--overwrite"),  type = "character", default = "f",
              help = "Overwrite previous output [default=%default]", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_bool <- function(x, flag) {
  if (tolower(x) %in% c("t", "true"))  return(TRUE)
  if (tolower(x) %in% c("f", "false")) return(FALSE)
  stop(paste0("--", flag, " must be t/f (got '", x, "')"), call. = FALSE)
}
SINGLE_END <- parse_bool(opt$single_end, "single_end")
OVERWRITE  <- parse_bool(opt$overwrite,  "overwrite")

if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  log_error("--input_dir and --output_dir are required arguments.")
  quit(status = 1)
}

INPUT_DIR  <- opt$input_dir
OUTPUT_DIR <- opt$output_dir
NSLOTS     <- opt$nslots
PATTERN_R1 <- opt$r1_pattern
PATTERN_R2 <- opt$r2_pattern

###############################################################################
### 3. Validate input directory and dependencies
###############################################################################

if (!dir.exists(INPUT_DIR)) {
  log_error(paste("Input directory does not exist:", INPUT_DIR))
  quit(status = 1)
}

registerDoParallel(cores = NSLOTS)

###############################################################################
### 4. Prepare output directories (output/, logs/, stats/)
###############################################################################

if (dir.exists(OUTPUT_DIR)) {
  if (!OVERWRITE) {
    log_error(paste("Output directory already exists:", OUTPUT_DIR,
                    "- use --overwrite t to overwrite"))
    quit(status = 1)
  }
  log_warn(paste("Overwriting existing directory:", OUTPUT_DIR))
  unlink(OUTPUT_DIR, recursive = TRUE)
}

RESULTS_DIR <- file.path(OUTPUT_DIR, "output")
LOGS_DIR    <- file.path(OUTPUT_DIR, "logs")
STATS_DIR   <- file.path(OUTPUT_DIR, "stats")
for (d in c(RESULTS_DIR, LOGS_DIR, STATS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_OUT   <- file.path(LOGS_DIR,  "1.2-quality-check.log")
STATS_OUT <- file.path(STATS_DIR, "1.2-quality-check-stats.tsv")

###############################################################################
### 5. Find input files
###############################################################################

rawR1 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R1, full.names = TRUE))
if (length(rawR1) == 0) {
  log_error(paste("No files found matching pattern:", PATTERN_R1))
  quit(status = 1)
}

if (!SINGLE_END) {
  rawR2 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R2, full.names = TRUE))
  if (length(rawR2) == 0) {
    log_error(paste("No R2 files found matching pattern:", PATTERN_R2))
    quit(status = 1)
  }
  if (length(rawR1) != length(rawR2)) {
    log_error(paste("Mismatch in number of files. R1:", length(rawR1),
                    "R2:", length(rawR2)))
    quit(status = 1)
  }
  log_msg(paste("Found", length(rawR1), "R1 files and", length(rawR2), "R2 files"))
} else {
  log_msg(paste("Found", length(rawR1), "single-end files"))
}

extract_sample_name <- function(filepath, pattern) {
  basename(filepath) %>% sub(pattern = pattern, replacement = "", fixed = FALSE)
}

###############################################################################
### 6. Per-sample read counts (R1 / single-end)
###############################################################################

count_seqs <- function(p, pattern) {
  n <- readFastq(p) %>% sread() %>% as.character() %>% length()
  data.frame(sample = extract_sample_name(p, pattern), nseq = n)
}

seq_counts_df <- foreach(i = rawR1, .combine = rbind) %dopar% {
  count_seqs(i, PATTERN_R1)
}

###############################################################################
### 7. R1 mean quality vs read count
###############################################################################

x_r1 <- qa(dirPath = INPUT_DIR, pattern = PATTERN_R1, sample = TRUE, n = 5000)
qa_means <- x_r1[["perCycle"]][["quality"]] %>%
            group_by(lane) %>%
            summarize(mean_q = sum(Score * Count) / sum(Count), .groups = "drop")
qa_means$lane <- sapply(qa_means$lane, function(x) extract_sample_name(x, PATTERN_R1))
qa_means2counts <- left_join(qa_means, seq_counts_df, by = c("lane" = "sample"))

text_size <- 2
p_r1 <- ggplot(qa_means2counts, aes(x = mean_q, y = nseq)) +
        geom_point() +
        scale_y_log10() +
        ylab("Read counts (log)") +
        xlab("Mean quality score (R1)") +
        geom_text(aes(label = as.character(lane)), hjust = 0.5, vjust = -1, size = text_size)

ggsave(p_r1, filename = file.path(RESULTS_DIR, "r1_mean_q_vs_nseq.png"),
       device = "png", width = 5, height = 4, dpi = 300)

# R1 mean quality kept for the stats table
r1_mean_q <- qa_means %>% rename(sample = lane, mean_q_r1 = mean_q)

###############################################################################
### 8. R2 mean quality vs read count (paired-end only)
###############################################################################

r2_mean_q <- NULL
if (!SINGLE_END) {
  x_r2 <- qa(dirPath = INPUT_DIR, pattern = PATTERN_R2, sample = TRUE, n = 5000)
  qa_means_r2 <- x_r2[["perCycle"]][["quality"]] %>%
                 group_by(lane) %>%
                 summarize(mean_q = sum(Score * Count) / sum(Count), .groups = "drop")
  qa_means_r2$lane <- sapply(qa_means_r2$lane, function(x) extract_sample_name(x, PATTERN_R2))
  qa_means2counts_r2 <- left_join(qa_means_r2, seq_counts_df, by = c("lane" = "sample"))

  p_r2 <- ggplot(qa_means2counts_r2, aes(x = mean_q, y = nseq)) +
          geom_point() +
          scale_y_log10() +
          ylab("Read counts (log)") +
          xlab("Mean quality score (R2)") +
          geom_text(aes(label = as.character(lane)), hjust = 0.5, vjust = -1, size = text_size)

  ggsave(p_r2, filename = file.path(RESULTS_DIR, "r2_mean_q_vs_nseq.png"),
         device = "png", width = 5, height = 4, dpi = 300)

  r2_mean_q <- qa_means_r2 %>% rename(sample = lane, mean_q_r2 = mean_q)
}

###############################################################################
### 9. Read-count histograms
###############################################################################

samples_hist_p <- ggplot(qa_means2counts, aes(nseq)) +
                  geom_histogram(bins = 40) + ylab("Num of samples")
samples_hist_log_p <- ggplot(qa_means2counts, aes(log(nseq))) +
                      geom_histogram(bins = 40) + ylab("Num of samples")

ggsave(samples_hist_p, filename = file.path(RESULTS_DIR, "samples_hist.png"),
       device = "png", width = 5, height = 4, dpi = 300)
ggsave(samples_hist_log_p, filename = file.path(RESULTS_DIR, "samples_hist_log.png"),
       device = "png", width = 5, height = 4, dpi = 300)

###############################################################################
### 10. PhiX contamination estimate
###############################################################################

count_phix_seqs <- function(p) {
  fastq <- readFastq(p) %>% sread() %>% as.character()
  fastq_nphix <- isPhiX(seqs = fastq, wordSize = 16, minMatches = 2)
  data.frame(file = p, n = length(fastq), nphix = sum(fastq_nphix))
}

phix_counts_df <- foreach(i = rawR1, .combine = rbind) %dopar% {
  count_phix_seqs(i)
}

X <- data.frame(
  sample = sapply(phix_counts_df$file, function(x) extract_sample_name(x, PATTERN_R1)),
  perc   = 100 * phix_counts_df$nphix / phix_counts_df$n,
  nphix  = phix_counts_df$nphix
)

X_plot <- X
undet <- grep(pattern = "Undetermined", X_plot$sample)
if (length(undet) > 0) X_plot <- X_plot[-undet, ]
X_plot$color <- ifelse(X_plot$perc > 0.005, "indianred", "gray50")

perc_phix_barplot <- ggplot(X_plot, aes(x = sample, y = perc, fill = color)) +
                     geom_bar(stat = "identity") +
                     theme_bw() +
                     scale_fill_manual(values = c("gray50" = "gray50", "indianred" = "indianred"),
                                       labels = c("gray50" = "Phix <= 0.005%", "indianred" = "Phix > 0.005%"),
                                       name = "") +
                     theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                           legend.position = "top")

ggsave(perc_phix_barplot, filename = file.path(RESULTS_DIR, "samples_perc_phix_barplot.png"),
       device = "png", width = 10, height = 4, dpi = 300)

###############################################################################
### 11. Stats table (samples as rows, statistics as columns)
###############################################################################

stats_tbl <- seq_counts_df %>%
             left_join(r1_mean_q, by = "sample")
if (!SINGLE_END && !is.null(r2_mean_q)) {
  stats_tbl <- stats_tbl %>% left_join(r2_mean_q, by = "sample")
}
stats_tbl <- stats_tbl %>%
             left_join(X %>% transmute(sample = sample, phix_pct = perc), by = "sample")

write_tsv(stats_tbl, STATS_OUT)

###############################################################################
### 12. Log and summary
###############################################################################

log_msg(paste("Processed", length(rawR1), "samples"))
log_msg("\033[0;32m1.2-quality-check.R completed successfully\033[0m")

generated <- c("r1_mean_q_vs_nseq.png",
               if (!SINGLE_END) "r2_mean_q_vs_nseq.png",
               "samples_hist.png", "samples_hist_log.png",
               "samples_perc_phix_barplot.png")

log_lines <- build_log(
  SCRIPT_NAME, SCRIPT_DESC, sample_name = "all samples",
  inputs  = c(paste("Input directory:", INPUT_DIR),
              paste("R1 pattern:", PATTERN_R1),
              if (!SINGLE_END) paste("R2 pattern:", PATTERN_R2)),
  params  = c(paste("Threads:", NSLOTS),
              paste("Read type:", if (SINGLE_END) "single-end" else "paired-end")),
  outputs = c(file.path(RESULTS_DIR, generated), paste("Statistics:", STATS_OUT)),
  command = paste(c(SCRIPT_NAME, commandArgs(trailingOnly = TRUE)), collapse = " "),
  exit_status = 0,
  tool_log = paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(log_lines, LOG_OUT)
