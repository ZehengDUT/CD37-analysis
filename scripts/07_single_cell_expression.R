# Summarize candidate-gene expression in annotated AAA and AS cells.
source(file.path("R", "common.R"))
require_packages(c("Seurat", "ggplot2"))
genes <- read_gene_list(project_path("config", "candidate_genes.txt"))
inputs <- c(aaa = project_path("data", "restricted", "aaa_annotated.rds"),
            as = project_path("data", "processed", "as_annotated.rds"))
for (cohort in names(inputs)) {
  if (!file.exists(inputs[[cohort]])) { message(cohort, ": annotated object missing; skip"); next }
  obj <- readRDS(inputs[[cohort]])
  if (!inherits(obj, "Seurat")) stop(cohort, " input is not a Seurat object")
  if (!"major_cell_type" %in% names(obj@meta.data)) stop(cohort, " requires major_cell_type metadata")
  Seurat::DefaultAssay(obj) <- "RNA"
  present <- intersect(genes, rownames(obj))
  if (!length(present)) stop(cohort, ": no candidate genes in RNA assay")
  labels <- as.character(obj@meta.data$major_cell_type)
  expression <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")[present, , drop = FALSE]
  groups <- split(seq_along(labels), labels)
  summary <- do.call(rbind, lapply(names(groups), function(type) {
    cells <- groups[[type]]
    data.frame(cell_type = type, gene = present, cells = length(cells),
               mean_expression = Matrix::rowMeans(expression[, cells, drop = FALSE]),
               percent_expressing = Matrix::rowMeans(expression[, cells, drop = FALSE] > 0) * 100)
  }))
  write_table(summary, project_path("results", "single_cell", paste0(cohort, "_gene_summary.tsv")))
  plot <- Seurat::DotPlot(obj, features = present, group.by = "major_cell_type") +
    ggplot2::coord_flip() + ggplot2::theme_bw()
  ggplot2::ggsave(project_path("results", "single_cell", paste0(cohort, "_candidate_dotplot.pdf")),
                  plot = plot, width = 11, height = 7)
}
