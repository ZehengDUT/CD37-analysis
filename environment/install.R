# Run with R 4.4.2. Package versions should be recorded with 99_session_info.R.
cran <- c("WGCNA", "glmnet", "caret", "e1071", "kernlab", "neuralnet", "pROC", "ggplot2", "readxl")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc <- c("limma", "sva", "clusterProfiler", "org.Hs.eg.db", "GEOquery")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
# CellChat 1.6.1, Seurat 4.3.0.1, and harmony 1.2.0 were present in the
# original analysis environment. Follow their upstream installation guides
# for matching versions. This script does not install an unspecified version.
