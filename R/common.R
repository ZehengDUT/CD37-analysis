repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file)) return(normalizePath(file.path(dirname(file[1]), ".."), mustWork = TRUE))
  normalizePath(".", mustWork = TRUE)
}

project_path <- function(...) file.path(repo_root(), ...)

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_expression <- function(path) {
  if (!file.exists(path)) stop("Missing expression matrix: ", path, call. = FALSE)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 5 || !names(x)[1] %in% c("gene", "id", "Gene", "symbol"))
    stop("Expected first column 'gene' or 'id' and at least four samples: ", path, call. = FALSE)
  genes <- trimws(as.character(x[[1]]))
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes))
    stop("Gene symbols must be nonempty and unique in ", path, call. = FALSE)
  m <- as.matrix(x[-1]); storage.mode(m) <- "double"; rownames(m) <- genes
  if (anyNA(m) || any(!is.finite(m))) stop("Expression matrix has missing or non-finite values: ", path, call. = FALSE)
  if (anyDuplicated(colnames(m))) stop("Duplicate sample names: ", path, call. = FALSE)
  m
}

sample_groups <- function(m) {
  s <- colnames(m)
  g <- ifelse(grepl("_con$", s), "Control", ifelse(grepl("_treat$", s), "Disease", NA))
  if (anyNA(g) || length(unique(g)) != 2L) stop("Sample names must end in _con or _treat and include both groups", call. = FALSE)
  factor(g, levels = c("Control", "Disease"))
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

read_gene_list <- function(path) {
  if (!file.exists(path)) stop("Missing gene list: ", path, call. = FALSE)
  unique(scan(path, what = character(), quiet = TRUE))
}
