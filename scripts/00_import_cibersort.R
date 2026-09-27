# Convert local CIBERSORTx workbooks into simple tables without publishing data.
# Save source workbooks as data/processed/{aaa,as}_cibersort_source.xlsx.
source(file.path("R", "common.R"))
require_packages("readxl")
for (cohort in c("aaa", "as")) {
  path <- project_path("data", "processed", paste0(cohort, "_cibersort_source.xlsx"))
  if (!file.exists(path)) { message(cohort, ": workbook missing; skip"); next }
  x <- as.data.frame(readxl::read_excel(path, skip = 10, .name_repair = "minimal"), check.names = FALSE)
  if (!all(c("Mixture", "P-value") %in% names(x))) stop(cohort, ": unexpected CIBERSORT workbook header")
  names(x)[names(x) == "Mixture"] <- "sample"
  names(x)[names(x) == "P-value"] <- "P.value"
  write_table(x, project_path("data", "processed", paste0(cohort, "_cibersort.tsv")))
  message(cohort, ": converted ", nrow(x), " samples")
}
