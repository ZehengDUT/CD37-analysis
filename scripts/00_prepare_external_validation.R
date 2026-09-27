# Convert an author-prepared, gene-symbol GSE47472 matrix to the input format.
# Place its source file at data/processed/gse47472_source.tsv; this stays local.
source(file.path("R", "common.R"))
path <- project_path("data", "processed", "gse47472_source.tsv")
if (!file.exists(path)) stop("Missing GSE47472 gene-symbol source matrix: ", path)
x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
if (!names(x)[1] %in% c("geneNames", "gene", "id")) stop("First column must hold gene symbols")
names(x)[1] <- "gene"
ids <- names(x)[-1]
cases <- paste0("GSM", 1150689:1150702)
controls <- paste0("GSM", 1150703:1150710)
if (!setequal(ids, c(cases, controls))) stop("GSE47472 columns do not match the 14 AAA and 8 donor accessions")
names(x)[-1] <- paste0(ids, ifelse(ids %in% cases, "_treat", "_con"))
if (anyDuplicated(x$gene)) stop("Duplicate gene symbols; resolve probe aggregation before use")
write_table(x, project_path("data", "processed", "aaa_validation.tsv"))
message("Prepared 14 AAA and 8 donor GSE47472 samples")
