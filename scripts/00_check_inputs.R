source(file.path("R", "common.R"))
for (cohort in c("aaa", "as")) {
  lasso <- read_gene_list(project_path("config", paste0("final_", cohort, "_lasso_genes.txt")))
  svm <- read_gene_list(project_path("config", paste0("final_", cohort, "_svm_rfe_genes.txt")))
  reported <- read_gene_list(project_path("config", paste0("final_", cohort, "_ml_genes.txt")))
  if (!setequal(intersect(lasso, svm), reported)) stop(cohort, " final ML lists disagree")
}
aaa <- read_gene_list(project_path("config", "final_aaa_ml_genes.txt"))
as <- read_gene_list(project_path("config", "final_as_ml_genes.txt"))
core <- read_gene_list(project_path("config", "final_core_genes.txt"))
if (length(core) != 5L || !setequal(intersect(aaa, as), core)) stop("Final shared gene list is inconsistent")
message("Final ML lists agree on: ", paste(core, collapse = ", "))
for (cohort in c("aaa", "as")) {
  path <- project_path("data", "processed", paste0(cohort, "_expression.tsv"))
  if (!file.exists(path)) { message(cohort, ": missing ", path); next }
  x <- read_expression(path); g <- sample_groups(x)
  message(cohort, ": ", nrow(x), " genes; ", ncol(x), " samples; ", paste(names(table(g)), table(g), collapse = ", "))
}
for (path in c("data/processed/aaa_validation.tsv", "data/restricted/aaa_annotated.rds"))
  message(path, if (file.exists(project_path(path))) ": present" else ": absent")
