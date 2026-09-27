# Apply author-reviewed cell labels before communication analysis.
source(file.path("R", "common.R"))
require_packages("Seurat")
path <- project_path("data", "restricted", "aaa_clustered.rds")
if (!file.exists(path)) stop("Missing clustered Seurat object: ", path)
obj <- readRDS(path)
labels_path <- project_path("data", "restricted", "aaa_cell_annotations.tsv")
if (!file.exists(labels_path)) stop("Missing reviewed cell labels: ", labels_path)
labels <- read.delim(labels_path, stringsAsFactors = FALSE)
required <- c("cell_id", "major_cell_type", "macrophage_subtype")
if (!all(required %in% names(labels)) || anyDuplicated(labels$cell_id))
  stop("Annotation table needs unique cell_id and both cell-type fields")
idx <- match(colnames(obj), labels$cell_id)
if (anyNA(idx)) stop("Missing annotations for ", sum(is.na(idx)), " cells")
obj$major_cell_type <- labels$major_cell_type[idx]
obj$macrophage_subtype <- labels$macrophage_subtype[idx]
saveRDS(obj, project_path("data", "restricted", "aaa_annotated.rds"))
message("Saved annotated object for CellChat")
