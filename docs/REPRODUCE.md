# Reproducing the analyses

## 1. Obtain the raw data

Download the 10x Genomics single-cell raw reads from **GSA PRJCA065583**
(controlled access; application via the Data Access Committee on the accession
page). Place the Cell Ranger output in `data/derived/` relative to the repository
root. The processed Seurat object used by the scripts is
`data/derived/seurat_reanalysis_v2_3grp.rds` (17,336 genes x 24,896 cells).

Larger intermediate objects (SCENIC database, montage previews, per-panel
exports) are not distributed; they can be regenerated from the scripts.

## 2. Environment

See `docs/ENVIRONMENT.md`. Core: R 4.4.2, Seurat 5.5.1, GSVA 2.0.7,
CellChat, Python 3 with pandas/matplotlib/Pillow.

```r
setwd("<this repository>")
source("scripts/_ssgsea_reproduce_and_fix.R")   # final ssGSEA table
```

All absolute paths were rewritten to repository-relative paths at release. If a
script still needs a machine-specific location, set it once at the top of the
script; no other changes are required.

## 3. Order of execution

| Step | Script | Output |
|---|---|---|
| 1 | `doublet_qc.R`, `make_S1_sample_umap_qc.R` | QC, Fig. S1 |
| 2 | `pipeline_r.R` | clustering, markers, proportions, DEG, LR table (see warning below) |
| 3 | `_ssgsea_reproduce_and_fix.R` | `tables/ssGSEA_pathway_activity.csv` |
| 4 | `make_figS9_v2.R` | Supplementary Fig. S9 |
| 5 | `cellchat_lr.R`, `cellchat_diff.R` | communication inference |
| 6 | `nichenet_analysis.R` | ligand-priority ranking |
| 7 | `make_fig3C_pathway_score.R`, `make_fig4B_tc_prop.R`, `_fig4_deg.R`, `_fig4E_func_v2.R` | main figure panels |
| 8 | `fig5ij_redraw.py`, `build_wet_archive.py`, `integrate_wet_scRNA.R` | wet-lab quantification |

> **Warning.** `pipeline_r.R` contains a legacy rank-mean pathway-activity block
> that is NOT the published ssGSEA result. Run step 3 separately and use its
> output. See `docs/METHODS_NOTES.md`.

## 4. Known limitations

* One animal per group; no biological replication at the animal level.
* Absolute viral-load values are reported as log10 genome copies per g of tissue
  (tissue) or per mL (serum); see the manuscript Methods for the exact n per assay.
* `nichenet_analysis.R` implements a heuristic dictionary score, not the
  `nichenetr` model; its ranking depends on the explicit directional prior
  defined in the script.
