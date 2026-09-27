source(file.path("R", "common.R"))
dir.create(project_path("results"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), project_path("results", "sessionInfo.txt"))
packages <- c("limma", "sva", "WGCNA", "glmnet", "caret", "e1071", "kernlab",
              "neuralnet", "pROC", "readxl", "clusterProfiler", "org.Hs.eg.db",
              "Seurat", "SeuratObject", "harmony", "CellChat", "ComplexHeatmap")
versions <- data.frame(package = packages,
                       version = vapply(packages, function(p)
                         if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_,
                         character(1)))
write_table(versions, project_path("results", "package_versions.tsv"))
message("Saved results/sessionInfo.txt")
