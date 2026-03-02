args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[[i + 1]]
}

ref <- get_arg("--ref")   # grch37 | grch38
#lib <- get_arg("--lib")   # path to workflow-local R library

## Setting additional libpath, but hopefullly it will work without it
# dir.create(lib, recursive = TRUE, showWarnings = FALSE)
# .libPaths(c(lib, .libPaths()))

# message("[CONFIG] ref=", ref)
# message("[CONFIG] lib=", lib)
# message("[R] .libPaths()=")
# print(.libPaths())

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

#options(repos = BiocManager::repositories())

# core packages
pkgs_bioc <- c(
  "MutationalPatterns",
  "VariantAnnotation",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "GenomeInfoDb",
  "Biostrings"
)

# ref-specific BSgenome
if (ref == "grch37") {
  pkgs_bioc <- c(pkgs_bioc, "BSgenome.Hsapiens.UCSC.hg19")
} else if (ref == "grch38") {
  pkgs_bioc <- c(pkgs_bioc, "BSgenome.Hsapiens.UCSC.hg38")
} else {
  stop("ref must be grch37 or grch38")
}

message("[INSTALL] Bioconductor packages:")
print(pkgs_bioc)


### Temporary workaround for bioc repo -------------------------------
# critical line: this controls BioCsoft/BioCann/... base URL
options(BioC_mirror = "https://ftp.gwdg.de/pub/misc/bioconductor")

# now build the full consistent repo set from BioC_mirror + CRAN
options(repos = BiocManager::repositories())

# optional: just silences the “getOption('repos') replaces …” message
options(BiocManager.check_repositories = FALSE)

print(getOption("repos"))

# ---------------------------------------------------------------------

BiocManager::install(pkgs_bioc, ask = FALSE, update = FALSE)

message("[CHECK] MutationalPatterns version: ", as.character(packageVersion("MutationalPatterns")))
if (ref == "grch37") message("[CHECK] hg19 version: ", as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg19")))
if (ref == "grch38") message("[CHECK] hg38 version: ", as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38")))

message("[DONE] prereqs installed into: ", .libPaths())
