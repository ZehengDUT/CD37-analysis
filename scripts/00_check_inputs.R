source(file.path("R", "common.R"))
for (cohort in c("aaa", "as")) {
  path <- project_path("data", "processed", paste0(cohort, "_expression.tsv"))
  if (!file.exists(path)) { message(cohort, ": missing ", path); next }
  x <- read_expression(path); g <- sample_groups(x)
  message(cohort, ": ", nrow(x), " genes; ", ncol(x), " samples; ", paste(names(table(g)), table(g), collapse = ", "))
}
for (path in c("data/processed/aaa_validation.tsv", "data/restricted/aaa_annotated.rds"))
  message(path, if (file.exists(project_path(path))) ": present" else ": absent")
