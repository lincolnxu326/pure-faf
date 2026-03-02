# -----------------------------------------------------
# Filter DeepSomatic intersected INDELs using muttype (MutationalPatterns)
# Rules:
#   - If VAF > 5%: DeepSomatic == PASS, Mutect2 == PASS
#   - If VAF > 10%: Mutect2 == PASS
#   - If Mutect2 == clustered_events: VAF > 30%
#   - If muttype ends with bp_insertion: drop unless (Mutect2 PASS and DeepSomatic PASS)
# -----------------------------------------------------

log <- file(snakemake@log[[1]], open = "at")
sink(log, append = TRUE)
sink(log, type = "message", append = TRUE)
message("-----------------------------------------------")
message("-----------------------------------------------")
message(date())

options(warn = 1)

suppressPackageStartupMessages({
  library(MutationalPatterns)
  library(VariantAnnotation)
  library(SummarizedExperiment)  # rowRanges generic lives here
  library(tibble)
  library(readr)
})

# ---- Compatibility patch: samples(CollapsedVCF) ----
# Some code (including some MutationalPatterns versions) calls samples(vcf) where vcf is CollapsedVCF.
# Define an S4 method so it resolves to samples(header(vcf)).
if (methods::isGeneric("samples") &&
    !methods::hasMethod("samples", signature(object = "CollapsedVCF"))) {
  methods::setMethod(
    "samples",
    signature(object = "CollapsedVCF"),
    function(object) VariantAnnotation::samples(VariantAnnotation::header(object))
  )
}

# ---- Reference genome ----
if (snakemake@params[["ref_genome"]] == "grch38") {
  suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
  })
  ref_genome <- BSgenome.Hsapiens.UCSC.hg38
} else if (snakemake@params[["ref_genome"]] == "grch37") {
  suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg19)
  })
  ref_genome <- BSgenome.Hsapiens.UCSC.hg19
} else {
  stop("ref_genome must be 'grch37' or 'grch38'")
}

# ---- Inputs ----
in_vcf  <- snakemake@input[["v"]]
out_vcf <- snakemake@output[["v"]]
out_rep <- snakemake@output[["report"]]
tum     <- snakemake@params[["tum"]]

message("[CONFIG] in=", in_vcf)
message("[CONFIG] out=", out_vcf)
message("[CONFIG] report=", out_rep)
message("[CONFIG] tum=", tum)
message("[CONFIG] ref_genome=", snakemake@params[["ref_genome"]])

# ---- Read VCF (keep both samples to preserve downstream concat compatibility) ----
vcf <- VariantAnnotation::readVcf(in_vcf)
n_in <- length(vcf)
message("[INPUT] n_records=", n_in)

# sample names are stored in the header
samps <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
tum_idx <- match(tum, samps)
if (is.na(tum_idx)) stop("Tumour sample not found in VCF header samples: ", tum)

# ---- Extract Mutect2 FILTER (CharacterList -> single string) ----
filter_list <- VariantAnnotation::fixed(vcf)$FILTER
m2_filter <- vapply(filter_list, function(x) {
  if (length(x) == 0 || all(is.na(x))) "."
  else paste(as.character(x), collapse = ";")
}, character(1))

# ---- Extract DeepSomatic FILTER annotation ----
if (!("DS_FILTER" %in% colnames(VariantAnnotation::info(vcf)))) {
  stop("INFO/DS_FILTER missing in input VCF.")
}
ds_raw <- VariantAnnotation::info(vcf)$DS_FILTER
ds_filter <- vapply(ds_raw, function(x) {
  if (length(x) == 0 || all(is.na(x)) || x == "") "."
  else as.character(x)[1]
}, character(1))


# ---- Extract tumour VAF from FORMAT/AF ----
g <- VariantAnnotation::geno(vcf)

vaf <- rep(NA_real_, length(vcf))

af_obj <- g$AF
if (!is.null(af_obj)) {
  if (is.matrix(af_obj)) {
    vaf <- suppressWarnings(as.numeric(af_obj[, tum_idx, drop = TRUE]))
  } else if (is.array(af_obj) && length(dim(af_obj)) == 3) {
    # sometimes stored as [variant, sample, allele]; take first ALT
    vaf <- suppressWarnings(as.numeric(af_obj[, tum_idx, 1, drop = TRUE]))
  } else {
    stop("Unexpected FORMAT/AF shape; expected matrix or 3D array.")
  }
} else {
  stop("FORMAT/AF missing; cannot extract tumour VAF here.")
}

# ---- Variant key for muttype mapping ----
rr <- SummarizedExperiment::rowRanges(vcf)
ref1 <- as.character(VariantAnnotation::ref(vcf))
alt1 <- vapply(VariantAnnotation::alt(vcf), function(x) as.character(x)[1], character(1))
key  <- paste0(as.character(GenomicRanges::seqnames(rr)), ":", BiocGenerics::start(rr), "_", ref1, "/", alt1)

# ---- MutationalPatterns indel context ----
# Use sample_names=tum so MutationalPatterns focuses on tumour genotype where relevant.
grl <- read_vcfs_as_granges(in_vcf, sample_names = tum, genome = ref_genome, type = "all")
indel_ctx_list <- get_indel_context(grl, ref_genome)
indel_ctx <- indel_ctx_list[[1]]

muttype_vec <- S4Vectors::mcols(indel_ctx)$muttype
names(muttype_vec) <- names(indel_ctx)

muttype <- unname(muttype_vec[key])
message("[muttype] mapped=", sum(!is.na(muttype)), " / ", length(muttype))

# ---- Apply rules ----
is_clustered <- grepl("clustered_events", m2_filter)
is_pon <- grepl("panel_of_normals", m2_filter)
is_pass <- (m2_filter == "PASS" | m2_filter == ".")

m2_is_passlike <- (is_pass | is_pon) & !is_clustered
m2_is_cluster  <- is_clustered
ds_is_pass <- (ds_filter == "PASS")

# Treat PASS and PON the same, but clustered overrides both
keep_base <-
  (m2_is_cluster  & (vaf > 0.30)) |
  (m2_is_passlike & (vaf > 0.10)) |
  (m2_is_passlike & ds_is_pass & (vaf > 0.05))

is_bp_insertion <- (!is.na(muttype)) & grepl("bp_insertion$", muttype)
keep <- keep_base & (!is_bp_insertion | (m2_is_passlike & ds_is_pass))
keep[is.na(keep)] <- FALSE

message("[FILTER] keep=", sum(keep), " drop=", sum(!keep))

# ---- Report ----
rep <- tibble(
  CHROM = as.character(GenomicRanges::seqnames(rr)),
  POS   = BiocGenerics::start(rr),
  REF   = ref1,
  ALT   = alt1,
  M2_FILTER = m2_filter,
  DS_FILTER = ds_filter,
  VAF = vaf,
  muttype = muttype,
  keep = keep
)
readr::write_tsv(rep, out_rep)

# ---- Write filtered VCF (bgzip + CSI) ----
vcf_f <- vcf[keep]
tmp_vcf <- tempfile(fileext = ".vcf")
VariantAnnotation::writeVcf(vcf_f, tmp_vcf)

cmd_bgzip <- sprintf("bgzip -f -c %s > %s", shQuote(tmp_vcf), shQuote(out_vcf))
cmd_index <- sprintf("bcftools index --csi -f %s", shQuote(out_vcf))

message("[WRITE] ", cmd_bgzip)
system(cmd_bgzip)

message("[INDEX] ", cmd_index)
system(cmd_index)

if (!file.exists(out_vcf) || file.info(out_vcf)$size == 0) {
  stop("Output VCF was not created or is empty: ", out_vcf)
}
message("[DONE] wrote ", out_vcf)
