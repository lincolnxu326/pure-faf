# MicroSEC_parameterised2.R
#
# MicroSEC: Microhomology-induced chimeric read-originating sequence
#           error cleaning pipeline for FFPE samples
#
# Snakemake script version:
#   - Uses snakemake@input, snakemake@output, snakemake@params, snakemake@threads.
#   - Assumes a single sample in sample_info TSV.
#   - Writes final .tsv.gz exactly to snakemake@output[["msec"]].
#   - Uses explicit tmp paths from snakemake@output (no tier guessing).
#   - No setwd(), everything is absolute paths.

# Expected Snakemake interface:
#   input:
#     sample_info:  TSV with 1 row:
#       V1: sample_name
#       V2: mutation_file
#       V3: bam_or_cram
#       V4: read_length
#       V5: adapter_1
#       V6+: optional organism/adapter_2/panel/reference_genome/simple_repeat_list
#   output:
#     msec:      final MicroSEC result .tsv.gz
#     slim_bam:  slim BAM (temp)
#     slim_bai:  slim BAM index (temp)
#     regions:   BED file with regions (temp)
#   params:
#     threshold_p:     numeric p-value threshold
#     mem_per_thread:  memory per thread in MB (e.g., 2000)
#     tool_log:        (optional) path for MicroSEC-specific log (will be appended)
#   threads:
#     integer, passed to samtools and MicroSEC

# ----------------- Logging -----------------

if ("tool_log" %in% names(snakemake@params)) {
  log_path <- snakemake@params[["tool_log"]]
} else {
  log_path <- snakemake@log[[1]]
}

dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
log <- file(log_path, open = "at")
sink(log, append = TRUE)
sink(log, type = "message", append = TRUE)
message("-----------------------------------------------")
message("-----------------------------------------------")
message(date())

close_all <- function(status = 0L) {
  message(date())
  message("-----------------------------------------------")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log), silent = TRUE)
  quit(save = "no", status = status)
}

##############################################
.libPaths(".Rlib")
suppressPackageStartupMessages({
  library(MicroSEC)
  library(dplyr)
  library(readr)
  library(stringr)
  library(Rsamtools)
  library(BiocGenerics)
  library(Biostrings)
  library(tools)
})

# ----------------- Snakemake inputs/outputs/params -----------------

# final output path (used as-is, absolute)
out_path <- normalizePath(snakemake@output[["msec"]], mustWork = FALSE)

# sample_info TSV
sample_list <- normalizePath(snakemake@input[["sample_info"]], mustWork = TRUE)

# params
progress_bar <- "Y"
threshold_p  <- as.numeric(snakemake@params[["threshold_p"]])
threads      <- as.integer(snakemake@threads)

# memory per thread (MB) -> samtools -m string
params_list <- snakemake@params
if ("mem_per_thread" %in% names(params_list)) {
  memory_per_thread <- paste0(as.character(params_list[["mem_per_thread"]]), "M")
} else {
  memory_per_thread <- "4000M"  # default 4GB if not set
}

# "wd" is only used to sanity-check the final output directory
wd <- dirname(out_path)
if (!dir.exists(wd)) {
  stop(sprintf("Working/output directory does not exist: %s", wd))
}

# ----------------- Helpers -----------------

run_cmd <- function(cmd) {
  status <- system(cmd)
  if (status != 0) stop(sprintf("Command failed (exit %s): %s", status, cmd))
}

qpath <- function(x) {
  # Normalise but do not require existence (for new files)
  shQuote(normalizePath(x, mustWork = FALSE))
}

# Header-only MicroSEC TSV (gz). This must match your microsec_annotate expected header.
write_microsec_header_only <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  hdr <- paste(
    c(
      "Sample","Mut_type","Chr","Pos","Ref","Alt","SimpleRepeat_TRF","Neighborhood_sequence",
      "percent_Alt","sum_AD","Chr_original","Chr_bam","read_length","total_read","soft_clipped_read",
      "flag_hairpin","pre_support_length","post_support_length","short_support_length","pre_farthest",
      "post_farthest","low_quality_base_rate_under_q18","low_quality_pre","low_quality_post",
      "distant_homology_rate","soft_clipped_rate","prob_filter_1","prob_filter_3_pre","prob_filter_3_post",
      "filter_1_mutation_intra_hairpin_loop","filter_2_hairpin_structure",
      "filter_3_microhomology_induced_mutation","filter_4_highly_homologous_region",
      "filter_5_soft_clipped_reads","filter_6_simple_repeat","filter_7_mutation_at_homopolymer",
      "filter_8_low_quality","msec_filter_123","msec_filter_1234","msec_filter_all","comment"
    ),
    collapse = "\t"
  )
  con <- gzfile(path, open = "wt")
  writeLines(hdr, con = con)
  close(con)
}

# "empty" meaning: no data lines beyond header (or file missing/0 bytes)
is_empty_mut_file <- function(f) {
  if (!file.exists(f)) return(TRUE)
  sz <- file.info(f)$size
  if (is.na(sz) || sz == 0L) return(TRUE)
  x <- readLines(f, n = 50L, warn = FALSE)
  x <- x[nzchar(trimws(x))]
  length(x) < 2L
}

# ----------------- Read and unpack sample_info (single sample) -----------------

sample_info <- read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
if (nrow(sample_info) != 1) {
  stop(sprintf("This script expects exactly one sample in sample_info, but found %d.", nrow(sample_info)))
}

sample_name   <- sample_info[1, 1]
mutation_file <- sample_info[1, 2]
bam_file_in   <- sample_info[1, 3]
read_length   <- as.integer(sample_info[1, 4])
adapter_1     <- sample_info[1, 5]

# Optional columns handling (same logic as original)
if (sample_info[1, 6] %in% c("Human","Mouse","hg19","hg38","mm10")) {
  adapter_2 <- adapter_1
  organism  <- sample_info[1, 6]
  panel     <- sample_info[1, 7]
  if (ncol(sample_info) == 8)  reference_genome <- sample_info[1, 8]
  if (ncol(sample_info) == 9){ reference_genome <- sample_info[1, 8]; simple_repeat_list <- sample_info[1, 9] }
} else {
  adapter_2 <- sample_info[1, 6]
  organism  <- sample_info[1, 7]
  panel     <- sample_info[1, 8]
  if (ncol(sample_info) == 9)  reference_genome <- sample_info[1, 9]
  if (ncol(sample_info) == 10){ reference_genome <- sample_info[1, 9]; simple_repeat_list <- sample_info[1,10] }
}
if (is.na(adapter_2) || adapter_2 == "") adapter_2 <- adapter_1

# Normalise paths that are used by external tools
mutation_file <- normalizePath(mutation_file, mustWork = TRUE)
bam_file_in   <- normalizePath(bam_file_in, mustWork = TRUE)
if (exists("reference_genome")) {
  reference_genome <- normalizePath(reference_genome, mustWork = TRUE)
}

# ----------------- Paths for intermediates (from Snakemake outputs) -----------------

bam_file_slim     <- normalizePath(snakemake@output[["slim_bam"]], mustWork = FALSE)
bam_file_slim_bai <- normalizePath(snakemake@output[["slim_bai"]], mustWork = FALSE)
bed_path          <- normalizePath(snakemake@output[["regions"]], mustWork = FALSE)

tmp_dir <- dirname(bam_file_slim)
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------- EARLY EXIT: empty mutation_info -----------------

if (is_empty_mut_file(mutation_file)) {
  message(sprintf("[%s] mutation_info is empty/header-only. Skipping MicroSEC and emitting header-only TSV + placeholder outputs.", sample_name))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(bed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(bam_file_slim), recursive = TRUE, showWarnings = FALSE)

  file.create(bed_path)
  file.create(bam_file_slim)
  file.create(bam_file_slim_bai)

  write_microsec_header_only(out_path)
  close_all(status = 0L)
}

# ----------------- Genome & chromosome naming -----------------

ref_genome <- fun_load_genome(organism)
chr_no     <- fun_load_chr_no(organism)

# ----------------- Load mutations (ORIGINAL MicroSEC loader) -----------------

df_mutation <- fun_load_mutation(mutation_file, sample_name, ref_genome, chr_no)

# If loader returns nothing, exit cleanly
if (nrow(df_mutation) == 0) {
  message(sprintf("[%s] fun_load_mutation returned 0 rows. Emitting header-only result.", sample_name))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(bed_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(bam_file_slim), recursive = TRUE, showWarnings = FALSE)

  file.create(bed_path)
  file.create(bam_file_slim)
  file.create(bam_file_slim_bai)

  write_microsec_header_only(out_path)
  close_all(status = 0L)
}

# ----------------- CRITICAL: normalise Chr naming to match BSgenome -----------------
# Prevent: Error in ref_genome[[chrom]] : no such sequence
ref_has_chr <- any(grepl("^chr", ref_genome@user_seqnames))

df_mutation <- df_mutation %>%
  mutate(
    Chr = as.character(Chr),
    Chr = if (ref_has_chr) {
      ifelse(grepl("^chr", Chr), Chr, paste0("chr", Chr))
    } else {
      sub("^chr", "", Chr)
    }
  ) %>%
  filter(Chr %in% ref_genome@user_seqnames)

if (nrow(df_mutation) == 0) {
  message(sprintf("[%s] No mutations remain after contig normalisation to BSgenome. Emitting header-only result.", sample_name))
  write_microsec_header_only(out_path)
  close_all(status = 0L)
}

# ----------------- Harmonise mutation chromosomes to BAM/CRAM space (separate column) -----------------

header_raw <- tryCatch(
  system(sprintf("bash -lc 'samtools view -H %s'", qpath(bam_file_in)), intern = TRUE),
  error = function(e) stop(sprintf("Failed to read BAM/CRAM header for %s via samtools: %s", bam_file_in, conditionMessage(e)))
)

sq_lines <- header_raw[grepl("^@SQ", header_raw)]
if (length(sq_lines) == 0) stop(sprintf("No @SQ lines found in BAM/CRAM header for %s", bam_file_in))

bam_targets <- sub(".*SN:([^ \t]+).*", "\\1", sq_lines)
bam_has_chr <- any(grepl("^chr", bam_targets))

df_mutation <- df_mutation %>%
  mutate(
    Chr_bam = dplyr::case_when(
      bam_has_chr && !ref_has_chr ~ ifelse(grepl("^chr", Chr), Chr, paste0("chr", Chr)),
      !bam_has_chr && ref_has_chr ~ sub("^chr", "", Chr),
      TRUE ~ Chr
    )
  )

message(sprintf("[%s] BAM/CRAM contigs example: %s", sample_name, paste(head(bam_targets, 5), collapse = ", ")))
message(sprintf("[%s] Unique mutation chromosomes (ref space): %s", sample_name, paste(sort(unique(df_mutation$Chr)), collapse = ", ")))
message(sprintf("[%s] Unique mutation chromosomes (BAM space): %s", sample_name, paste(sort(unique(df_mutation$Chr_bam)), collapse = ", ")))

# Filter to BAM contigs BEFORE building BED (prevents empty BED/samtools weirdness)
df_mutation <- df_mutation %>% filter(Chr_bam %in% bam_targets)
if (nrow(df_mutation) == 0) {
  message(sprintf("[%s] No mutations are on BAM contigs after subsetting (pre-BED). Emitting header-only result.", sample_name))
  write_microsec_header_only(out_path)
  close_all(status = 0L)
}

# ----------------- Build BED around mutations (±200; merge if next within 400bp) -----------------

df_mutation_save <- df_mutation
download_region <- data.frame(matrix(rep(NA,3), nrow=1))[numeric(0),]
colnames(download_region) <- c("chrom","chromStart","chromEnd")

# Only chromosomes that exist both in BAM and in the mutations (BAM space)
chromosomes_bam <- intersect(bam_targets, unique(df_mutation_save$Chr_bam))

for (i in chromosomes_bam) {
  df_chr <- df_mutation_save[df_mutation_save$Chr_bam == i, ]
  if (nrow(df_chr) == 0) next
  start <- max(1, df_chr$Pos[1] - 200)
  for (k in seq_len(nrow(df_chr))) {
    if (k == nrow(df_chr) || df_chr$Pos[k+1] - df_chr$Pos[k] > 400) {
      download_region <- rbind(download_region, c(i, start - 1, df_chr$Pos[k] + 200))
      if (k < nrow(df_chr)) start <- max(1, df_chr$Pos[k+1] - 200)
    }
  }
}

dir.create(dirname(bed_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(download_region, bed_path, col_names = FALSE, progress = FALSE)

# ----------------- samtools pipeline (streamed) -----------------

mem_per_thread <- memory_per_thread
sort_prefix <- file.path(tmp_dir, "samtools_sort")
input_ext <- tolower(file_ext(bam_file_in))

if (input_ext == "bam") {
  cmd <- sprintf(
    paste(
      "bash -lc 'set -euo pipefail;",
      "export TMPDIR=%s;",
      "samtools view -@ %d -b -h --no-PG %s -L %s |",
      "samtools sort -@ %d -m %s -O BAM -T %s -o %s -'",
      sep=" "
    ),
    qpath(tmp_dir), threads, qpath(bam_file_in), qpath(bed_path),
    threads, mem_per_thread, qpath(sort_prefix), qpath(bam_file_slim)
  )
} else if (input_ext == "cram") {
  if (!exists("reference_genome")) stop("CRAM input requires reference_genome in sample info TSV.")
  cmd <- sprintf(
    paste(
      "bash -lc 'set -euo pipefail;",
      "export TMPDIR=%s;",
      "samtools view -@ %d -b -h --no-PG %s -T %s -L %s |",
      "samtools sort -@ %d -m %s -O BAM -T %s -o %s -'",
      sep=" "
    ),
    qpath(tmp_dir), threads, qpath(bam_file_in), qpath(reference_genome), qpath(bed_path),
    threads, mem_per_thread, qpath(sort_prefix), qpath(bam_file_slim)
  )
} else {
  stop(sprintf("Unsupported input extension: %s", input_ext))
}

message("Trimming and sorting reads into slim BAM (tmp/)...")
run_cmd(cmd)
run_cmd(sprintf("bash -lc 'set -euo pipefail; samtools index %s'", qpath(bam_file_slim)))

# ----------------- load slim BAM -----------------

df_bam <- fun_load_bam(bam_file_slim)

# ----------------- Filter mutations to slim BAM contigs (BAM space) -----------------

targets <- names(Rsamtools::scanBamHeader(bam_file_slim)[[1]]$targets)
missing_chr <- setdiff(unique(df_mutation$Chr_bam), targets)
if (length(missing_chr) > 0) {
  message(sprintf("[%s] Skipping contigs not in slim BAM (BAM space): %s", sample_name, paste(missing_chr, collapse = ", ")))
}
df_mutation <- dplyr::filter(df_mutation, Chr_bam %in% targets)

if (nrow(df_mutation) == 0) {
  message(sprintf("[%s] No mutations left on slim BAM contigs after subsetting; writing header-only result.", sample_name))
  write_microsec_header_only(out_path)
  close_all(status = 0L)
}

# ----------------- MicroSEC core (single sample) -----------------

result <- fun_read_check(
  df_mutation = df_mutation,
  df_bam = df_bam,
  ref_genome = ref_genome,
  sample_name = sample_name,
  read_length = read_length,
  adapter_1 = adapter_1,
  adapter_2 = adapter_2,
  short_homology_search_length = 4,
  min_homology_search = 40,
  progress_bar = progress_bar
)

msec            <- result[[1]]
homology_search <- result[[2]]
mut_depth       <- result[[3]]

msec <- fun_homology(
  msec,
  homology_search,
  min_homology_search = 40,
  ref_genome,
  chr_no,
  progress_bar = progress_bar
)
msec <- fun_summary(msec)
msec <- fun_analysis(
  msec,
  mut_depth,
  short_homology_search_length = 4,
  min_homology_search = 40,
  threshold_p = threshold_p,
  threshold_hairpin_ratio = 0.50,
  threshold_short_length = 0.75,
  threshold_distant_homology = 0.15,
  threshold_soft_clip_ratio = 0.50,
  threshold_low_quality_rate = 0.1,
  homopolymer_length = 15
)

# ----------------- Save final result -----------------

fun_save(msec, out_path)
message(sprintf("MicroSEC finished successfully. Output: %s", out_path))

close_all(status = 0L)