# AAA macrophage-centred secreted-signaling analysis with CellChat 1.6.1.
# Input stays in data/restricted/ and is never committed.
source(file.path("R", "common.R"))
require_packages(c("Seurat", "CellChat", "ComplexHeatmap"))
path <- project_path("data", "restricted", "aaa_annotated.rds")
if (!file.exists(path)) stop("Missing controlled-access annotated AAA Seurat object: ", path, call. = FALSE)
obj <- readRDS(path)
if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object")
required <- c("sample_id", "disease_status", "major_cell_type", "macrophage_subtype")
if (!all(required %in% names(obj@meta.data)))
  stop("Seurat metadata must include: ", paste(required, collapse = ", "), call. = FALSE)
obj <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data$disease_status == "AAA"])
if (!ncol(obj)) stop("No AAA cells found in disease_status metadata")
Seurat::DefaultAssay(obj) <- "RNA"
allowed_macrophages <- c("Inflammatory", "APOC1+ lipid-associated", "ISG+ SPP1+", "FOLR2+ resident-like")
major <- as.character(obj@meta.data$major_cell_type)
subtype <- as.character(obj@meta.data$macrophage_subtype)
label <- major
is_mac <- major == "Macrophage"
label[is_mac] <- subtype[is_mac]
invalid <- is.na(label) | !nzchar(label) | label %in% c("Doublet", "Unassigned", "Unknown") |
           (is_mac & !label %in% allowed_macrophages)
obj <- subset(obj, cells = colnames(obj)[!invalid])
obj$cellchat_label <- label[!invalid]
counts <- table(obj$cellchat_label)
keep <- names(counts)[counts >= 10L]
obj <- subset(obj, cells = colnames(obj)[obj$cellchat_label %in% keep])
if (!all(allowed_macrophages %in% obj$cellchat_label))
  stop("All four macrophage subtypes need at least 10 AAA cells after filtering")
if (!ncol(obj)) stop("No cells remain after group-size filter")
cellchat <- CellChat::createCellChat(object = obj, group.by = "cellchat_label", assay = "RNA")
data("CellChatDB.human", package = "CellChat", envir = environment())
database <- get("CellChatDB.human", envir = environment())
cellchat@DB <- CellChat::subsetDB(database, search = "Secreted Signaling")
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)
cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
out <- project_path("results", "cellchat")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(cellchat, file.path(out, "cellchat_aaa.rds"))
write_table(data.frame(cell_type = names(table(obj$cellchat_label)),
                       cells = as.vector(table(obj$cellchat_label))), file.path(out, "group_counts.tsv"))
all_interactions <- CellChat::subsetCommunication(cellchat)
write_table(all_interactions, file.path(out, "all_secreted_interactions.tsv"))
mac <- intersect(allowed_macrophages, levels(cellchat@idents))
outgoing <- CellChat::subsetCommunication(cellchat, sources.use = mac)
incoming <- CellChat::subsetCommunication(cellchat, targets.use = mac)
write_table(outgoing, file.path(out, "macrophage_outgoing.tsv"))
write_table(incoming, file.path(out, "macrophage_incoming.tsv"))
plot_pdf <- function(filename, plot_fun, width = 10, height = 8) {
  grDevices::pdf(file.path(out, filename), width = width, height = height, onefile = TRUE)
  tryCatch(plot_fun(), finally = grDevices::dev.off())
}
plot_pdf("global_interaction_count.pdf", function() {
  CellChat::netVisual_circle(cellchat@net$count, weight.scale = TRUE,
                             vertex.weight = as.vector(table(cellchat@idents)))
})
plot_pdf("global_interaction_strength.pdf", function() {
  CellChat::netVisual_circle(cellchat@net$weight, weight.scale = TRUE,
                             vertex.weight = as.vector(table(cellchat@idents)))
})
plot_pdf("macrophage_outgoing_bubble.pdf", function()
  print(CellChat::netVisual_bubble(cellchat, sources.use = mac, remove.isolate = TRUE)), 14, 8)
plot_pdf("macrophage_incoming_bubble.pdf", function()
  print(CellChat::netVisual_bubble(cellchat, targets.use = mac, remove.isolate = TRUE)), 14, 8)
focus <- c("SPP1", "MIF", "CCL", "CXCL", "TNF", "VEGF", "ANNEXIN", "GALECTIN")
available <- intersect(focus, cellchat@netP$pathways)
for (pathway in available) {
  safe <- gsub("[^A-Za-z0-9]", "_", pathway)
  write_table(CellChat::subsetCommunication(cellchat, signaling = pathway),
              file.path(out, paste0(safe, "_interactions.tsv")))
  plot_pdf(paste0(safe, "_circle.pdf"), function()
    CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "circle"))
  plot_pdf(paste0(safe, "_outgoing_heatmap.pdf"), function()
    ComplexHeatmap::draw(CellChat::netAnalysis_signalingRole_heatmap(cellchat, signaling = pathway,
                                                                      pattern = "outgoing")))
  plot_pdf(paste0(safe, "_incoming_heatmap.pdf"), function()
    ComplexHeatmap::draw(CellChat::netAnalysis_signalingRole_heatmap(cellchat, signaling = pathway,
                                                                      pattern = "incoming")))
}
writeLines(c("Requested pathways:", focus, "", "Detected pathways:", available),
           file.path(out, "pathway_status.txt"))
message("CellChat complete. Detected focus pathways: ", paste(available, collapse = ", "))
