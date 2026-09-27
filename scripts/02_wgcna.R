# Run WGCNA separately for AAA and AS from processed matrices.
source(file.path("R", "common.R"))
require_packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads(nThreads = 2)
for (cohort in c("aaa", "as")) {
  expr <- read_expression(project_path("data", "processed", paste0(cohort, "_expression.tsv")))
  group <- sample_groups(expr)
  variable <- apply(expr, 1, stats::var)
  keep <- order(variable, decreasing = TRUE)[seq_len(ceiling(length(variable) * 0.25))]
  dat <- t(expr[keep, , drop = FALSE])
  qc <- WGCNA::goodSamplesGenes(dat, verbose = 0)
  dat <- dat[qc$goodSamples, qc$goodGenes, drop = FALSE]
  trait <- as.integer(group[qc$goodSamples] == "Disease")
  sample_tree <- stats::hclust(stats::dist(dat), method = "average")
  # A height of 20,000 was reported. With log-scale expression it may exclude no sample.
  cluster <- WGCNA::cutreeStatic(sample_tree, cutHeight = 20000, minSize = 10)
  if (any(cluster == 0)) {
    dat <- dat[cluster != 0, , drop = FALSE]
    trait <- trait[cluster != 0]
  }
  if (nrow(dat) < 10L) stop(cohort, ": fewer than 10 samples after outlier exclusion")
  powers <- 1:30
  sft <- WGCNA::pickSoftThreshold(dat, powerVector = powers, verbose = 0)
  fit <- sft$fitIndices
  power <- unname(c(aaa = 15L, as = 7L)[cohort])
  if (is.na(power)) stop("No manuscript soft power for ", cohort)
  net <- WGCNA::blockwiseModules(dat, power = power, networkType = "unsigned",
                                 TOMType = "unsigned", minModuleSize = 60,
                                 mergeCutHeight = 0.25, maxBlockSize = ncol(dat),
                                 numericLabels = FALSE, saveTOMs = FALSE,
                                 verbose = 0)
  colors <- net$colors; names(colors) <- colnames(dat)
  mes <- net$MEs
  module_cor <- stats::cor(mes, trait, use = "pairwise.complete.obs")
  module_p <- WGCNA::corPvalueStudent(module_cor, nrow(dat))
  gs <- as.vector(stats::cor(dat, trait, use = "pairwise.complete.obs"))
  mm <- stats::cor(dat, mes, use = "pairwise.complete.obs")
  own_mm <- vapply(seq_along(colors), function(i) {
    col <- paste0("ME", colors[i]); if (col %in% colnames(mm)) mm[i, col] else NA_real_
  }, numeric(1))
  genes <- data.frame(gene = colnames(dat), module = unname(colors),
                      GS = gs, MM = own_mm, stringsAsFactors = FALSE)
  genes$hub_manuscript_threshold <- abs(genes$GS) > 0.6 & abs(genes$MM) > 0.8 & genes$module != "grey"
  prefix <- project_path("results", "wgcna", cohort)
  write_table(fit, paste0(prefix, "_soft_threshold.tsv"))
  write_table(data.frame(module = rownames(module_cor), correlation = as.vector(module_cor),
                         p_value = as.vector(module_p)), paste0(prefix, "_modules.tsv"))
  write_table(genes, paste0(prefix, "_genes.tsv"))
  observed_r2 <- fit$SFT.R.sq[match(power, fit$Power)]
  writeLines(c(paste0("soft_power=", power), paste0("observed_scale_free_R2=", observed_r2)),
             paste0(prefix, "_settings.txt"))
  if (is.finite(observed_r2) && observed_r2 < 0.90)
    warning(cohort, ": manuscript power ", power, " does not reach R^2 >= 0.90 in this input matrix")
  message(cohort, ": power ", power, ", ", sum(genes$hub_manuscript_threshold), " genes pass GS/MM cutoffs")
}
