# Input: processed, normalized gene × sample matrices with _con/_treat sample suffixes.
source(file.path("R", "common.R"))
require_packages("limma")
for (cohort in c("aaa", "as")) {
  expr <- read_expression(project_path("data", "processed", paste0(cohort, "_expression.tsv")))
  group <- sample_groups(expr)
  design <- stats::model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(expr, design))
  tab <- limma::topTable(fit, coef = "groupDisease", number = Inf, sort.by = "none")
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
  tab$selected <- tab$adj.P.Val < 0.05 & abs(tab$logFC) > 1
  write_table(tab, project_path("results", "de", paste0(cohort, "_all.tsv")))
  write_table(tab[tab$selected, ], project_path("results", "de", paste0(cohort, "_selected.tsv")))
  message(cohort, ": ", sum(tab$selected), " DE genes")
}
