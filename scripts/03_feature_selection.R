# Recalculate selection from the 18-gene DEG/WGCNA intersection.
# LASSO is binomial with alpha=1. SVM-RFE is a regression SVM evaluated by RMSE,
# matching the manuscript's stated selection metric (classification is assessed later).
source(file.path("R", "common.R"))
require_packages(c("glmnet", "caret", "kernlab"))
genes <- read_gene_list(project_path("config", "candidate_genes.txt"))
if (length(genes) != 18L) stop("Expected exactly 18 candidate genes")
selected <- list()
for (cohort in c("aaa", "as")) {
  expr <- read_expression(project_path("data", "processed", paste0(cohort, "_expression.tsv")))
  if (!all(genes %in% rownames(expr))) stop(cohort, " expression is missing: ", paste(setdiff(genes, rownames(expr)), collapse = ", "))
  x <- t(expr[genes, , drop = FALSE]); y <- as.integer(sample_groups(expr) == "Disease")
  lasso <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10, type.measure = "deviance")
  b <- as.matrix(stats::coef(lasso, s = "lambda.min"))
  lasso_genes <- setdiff(rownames(b)[as.vector(b != 0)], "(Intercept)")
  control <- caret::rfeControl(functions = caret::caretFuncs, method = "cv", number = 5,
                               returnResamp = "final", verbose = FALSE)
  svm <- caret::rfe(x = as.data.frame(x), y = y, sizes = seq_along(genes),
                    rfeControl = control, method = "svmRadial", metric = "RMSE",
                    trControl = caret::trainControl(method = "cv", number = 5),
                    tuneLength = 3)
  svm_genes <- as.character(caret::predictors(svm))
  both <- intersect(lasso_genes, svm_genes)
  out <- data.frame(gene = genes, LASSO = genes %in% lasso_genes,
                    SVM_RFE = genes %in% svm_genes, both = genes %in% both)
  write_table(out, project_path("results", "features", paste0(cohort, "_selection.tsv")))
  writeLines(both, project_path("results", "features", paste0(cohort, "_selected.txt")))
  saveRDS(list(lasso = lasso, svm = svm), project_path("results", "features", paste0(cohort, "_models.rds")))
  selected[[cohort]] <- both
  message(cohort, ": LASSO ", length(lasso_genes), ", SVM-RFE ", length(svm_genes), ", intersection ", length(both))
}
shared <- intersect(selected$aaa, selected$as)
writeLines(shared, project_path("results", "features", "shared_core_genes.txt"))
message("Shared across AAA and AS: ", paste(shared, collapse = ", "))
