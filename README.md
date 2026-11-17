# AAA–AS immune transcriptomics: reproducible analysis code

This repository provides R analysis scripts, documentation, configuration files, and final feature-gene lists for investigating CD37 and shared immune-related features in abdominal aortic aneurysm (AAA) and atherosclerosis (AS). The workflow covers differential expression, WGCNA, feature selection, diagnostic modelling, immune deconvolution, and AAA macrophage-centred CellChat analysis. Patient-level data and third-party commercial scripts are not redistributed.

## Reproducibility status

The code is runnable with the input formats described below. The corrected final AAA and AS machine-learning gene lists supplied by the authors are in `config/final_aaa_ml_genes.txt` and `config/final_as_ml_genes.txt`; their shared set is in `config/final_core_genes.txt`. `03_feature_selection.R` independently recalculates selection and reports any difference without overwriting the final lists. Exact numerical results should be checked against the final manuscript before archiving a release. Some AAA single-cell inputs are controlled-access under HRA011745 and must be obtained through the repository's access process.

## Public and controlled inputs

| Input | Source | Local destination | Format expected by code |
| --- | --- | --- | --- |
| AS bulk | [GSE100927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100927) | `data/processed/as_expression.tsv` | gene × sample, `gene` first column; sample names end `_con` or `_treat` |
| AAA bulk | [GSE183464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183464) | part of `data/processed/aaa_expression.tsv` | gene × sample; same suffix convention |
| AAA external validation | [GSE47472](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47472) | `data/processed/aaa_validation.tsv` | gene × sample; same suffix convention |
| AAA/controls single cell | [HRA011745](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011745) | `data/restricted/aaa_annotated.rds` | annotated Seurat object; metadata below |
| AS single cell | [GSE253903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253903) | local only | used for descriptive mapping; not required for the scripts below |

The processed expression matrices in the author's current folders can be used for a local audit: `AAA/normalize.txt` maps to `aaa_expression.tsv`, and `AS/normalize.txt` maps to `as_expression.tsv`. `GEO/GSE47472/GSE47472.txt` can be copied locally as `gse47472_source.tsv` and converted with `00_prepare_external_validation.R`. The existing CIBERSORTx workbooks can be copied locally as `{aaa,as}_cibersort_source.xlsx` and converted with `00_import_cibersort.R`. None of these data files are committed. Rebuilding the exact AAA and AS matrices from raw records requires the original normalization and batch-correction decisions, which cannot be inferred unambiguously from the supplied files.

The AAA CellChat script expects a Seurat object with a normalized `RNA` assay and metadata columns `sample_id`, `major_cell_type`, and `macrophage_subtype`. The latter must distinguish `Inflammatory`, `APOC1+ lipid-associated`, `ISG+ SPP1+`, and `FOLR2+ resident-like` for macrophages. Other cells use the broad labels in the manuscript. `data/restricted/aaa_annotated.rds` is never published.

## Run

Use R 4.4.2. Install packages listed in [environment/install.R](environment/install.R), then run from the repository root:

```text
Rscript scripts/00_check_inputs.R
Rscript scripts/00_prepare_external_validation.R
Rscript scripts/00_import_cibersort.R
Rscript scripts/01_differential_expression.R
Rscript scripts/02_wgcna.R
Rscript scripts/02_enrichment.R
Rscript scripts/03_feature_selection.R
Rscript scripts/04_ann_validation.R
Rscript scripts/05_immune_correlations.R
Rscript scripts/06_prepare_aaa_single_cell.R
Rscript scripts/06_apply_aaa_annotations.R
Rscript scripts/06_cellchat_aaa.R
Rscript scripts/07_single_cell_expression.R
Rscript scripts/99_session_info.R
```

Each script has a short input contract and writes under `results/`, which Git excludes. Run the scripts in order. ANN and immune correlations use the five genes in `config/final_core_genes.txt`. Immune correlations require a CIBERSORT or CIBERSORTx result table supplied separately. CellChat requires the controlled AAA Seurat object and author-verified macrophage annotations.

The package intentionally stops with a clear error when required data, metadata, or a resolved gene list is missing. It does not manufacture placeholder results.

No random seed is set by these scripts. Cross-validation folds and neural-network fits may therefore vary between runs; exact numerical replication is not guaranteed.

## Sources and attribution

The code uses R/Bioconductor packages cited in the manuscript. CIBERSORT source code and its signature matrix are not redistributed. Users should obtain and run that software under its own terms, or supply an equivalent result table and state which tool generated it. The CellChat implementation is based on the public package API, not copied source. Record exact package versions with `scripts/99_session_info.R`.

## Release and citation

The source repository is available at [ZehengDUT/CD37-analysis](https://github.com/ZehengDUT/CD37-analysis). Archived releases are preserved in Zenodo. Cite the version-specific DOI and release version corresponding to the archive used, together with the repository URL.

An archived release is a fixed snapshot: later changes to this repository do not update an existing Zenodo archive. Updated source files must be archived as a new version before the corresponding archive DOI can be cited. Check numerical results against the study's reported results before using a release for exact replication.

Neither this repository nor its code archives include the controlled-access single-cell sequencing data. Access to those data must be requested through [HRA011745](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011745).
