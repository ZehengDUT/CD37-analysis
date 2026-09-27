# Five-gene ANN with cohort-wise median scores and external validation.
source(file.path("R", "common.R"))
require_packages(c("neuralnet", "pROC"))
genes <- read_gene_list(project_path("config", "final_core_genes.txt"))
if (length(genes) != 5L)
  stop("ANN requires five genes in config/final_core_genes.txt", call. = FALSE)
train <- read_expression(project_path("data", "processed", "aaa_expression.tsv"))
test <- read_expression(project_path("data", "processed", "aaa_validation.tsv"))
if (!all(genes %in% rownames(train)) || !all(genes %in% rownames(test))) stop("Core genes absent from training or validation matrix")
de <- read.delim(project_path("results", "de", "aaa_all.tsv"), check.names = FALSE)
direction <- setNames(sign(de$logFC), de$gene)
if (anyNA(direction[genes]) || any(direction[genes] == 0)) stop("Core gene directions cannot be determined from DE results")
score <- function(expr) {
  x <- t(expr[genes, , drop = FALSE])
  medians <- apply(expr[genes, , drop = FALSE], 1, stats::median)
  out <- sapply(genes, function(g) if (direction[g] > 0) as.integer(x[, g] > medians[g]) else as.integer(x[, g] < medians[g]))
  as.data.frame(out, check.names = FALSE)
}
x_train <- score(train); x_test <- score(test)
x_train$Disease <- as.integer(sample_groups(train) == "Disease")
formula <- stats::reformulate(genes, response = "Disease")
model <- neuralnet::neuralnet(formula, data = x_train, hidden = 5, linear.output = FALSE,
                              act.fct = "logistic", stepmax = 1e6)
predict_probability <- function(x) as.vector(neuralnet::compute(model, x[, genes, drop = FALSE])$net.result)
evaluate <- function(expr, x, cohort) {
  actual <- as.integer(sample_groups(expr) == "Disease")
  probability <- predict_probability(x)
  roc <- pROC::roc(actual, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  pred <- data.frame(sample = colnames(expr), actual = actual, probability = probability,
                     predicted = as.integer(probability >= 0.5))
  write_table(pred, project_path("results", "ann", paste0(cohort, "_predictions.tsv")))
  data.frame(cohort = cohort, n = nrow(pred), AUC = as.numeric(pROC::auc(roc)),
             accuracy = mean(pred$predicted == pred$actual))
}
metrics <- rbind(evaluate(train, x_train, "aaa_training"), evaluate(test, x_test, "gse47472_validation"))
write_table(metrics, project_path("results", "ann", "metrics.tsv"))
saveRDS(model, project_path("results", "ann", "model.rds"))
print(metrics)
