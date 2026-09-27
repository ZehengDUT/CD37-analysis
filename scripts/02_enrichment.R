# GO and KEGG enrichment for the 18 shared DEG/WGCNA candidates.
source(file.path("R", "common.R"))
require_packages(c("clusterProfiler", "org.Hs.eg.db"))
genes <- read_gene_list(project_path("config", "candidate_genes.txt"))
mapped <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db::org.Hs.eg.db)
if (!nrow(mapped)) stop("No candidate genes mapped to Entrez IDs")
ids <- unique(mapped$ENTREZID)
go <- clusterProfiler::enrichGO(gene = ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05, readable = TRUE)
kegg <- clusterProfiler::enrichKEGG(gene = ids, organism = "hsa",
                                    pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write_table(mapped, project_path("results", "enrichment", "gene_mapping.tsv"))
write_table(as.data.frame(go), project_path("results", "enrichment", "go.tsv"))
write_table(as.data.frame(kegg), project_path("results", "enrichment", "kegg.tsv"))
