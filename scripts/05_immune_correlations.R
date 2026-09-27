# Input: CIBERSORT/CIBERSORTx output with sample, P.value, and immune fractions.
source(file.path("R", "common.R"))
for (cohort in c("aaa", "as")) {
  path <- project_path("data", "processed", paste0(cohort, "_cibersort.tsv"))
  if (!file.exists(path)) stop("Missing independent CIBERSORT result: ", path, call. = FALSE)
  fractions <- read.delim(path, check.names = FALSE)
  if (!all(c("sample", "P.value") %in% names(fractions))) stop("CIBERSORT table requires sample and P.value")
  fractions <- fractions[is.finite(fractions$P.value) & fractions$P.value < 0.05, , drop = FALSE]
  expr <- read_expression(project_path("data", "processed", paste0(cohort, "_expression.tsv")))
  fractions <- fractions[match(colnames(expr), fractions$sample), , drop = FALSE]
  fractions <- fractions[!is.na(fractions$sample), , drop = FALSE]
  if (nrow(fractions) < 5) stop(cohort, ": fewer than five samples pass P<0.05")
  immune <- setdiff(names(fractions), c("sample", "P.value", "Correlation", "RMSE"))
  immune <- immune[vapply(fractions[immune], is.numeric, logical(1))]
  genes <- read_gene_list(project_path("results", "features", "shared_core_genes.txt"))
  if (length(genes) != 5L) stop("Immune correlation step requires five verified core genes")
  genes <- intersect(genes, rownames(expr))
  grid <- expand.grid(gene = genes, cell_type = immune, stringsAsFactors = FALSE)
  grid$rho <- NA_real_; grid$p_value <- NA_real_; grid$n <- nrow(fractions)
  for (i in seq_len(nrow(grid))) {
    z <- suppressWarnings(stats::cor.test(as.numeric(expr[grid$gene[i], fractions$sample]),
                                           fractions[[grid$cell_type[i]]], method = "spearman", exact = FALSE))
    grid$rho[i] <- unname(z$estimate); grid$p_value[i] <- z$p.value
  }
  grid$FDR <- stats::p.adjust(grid$p_value, method = "BH")
  write_table(grid, project_path("results", "immune", paste0(cohort, "_gene_correlations.tsv")))
}
