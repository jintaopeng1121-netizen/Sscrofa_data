# Sscrofa_data

Analysis code and derived tables for:

> **A multi-component herbal formulation limits PRRSV replication and reduces PD-L1-high neutrophil abundance in porcine lungs**
> Jintao Peng, Xu Wang, Wei Zhou, Xin Fan, Dishi Chen, Shijie Li, Pinpin Chen, Hongxia Ding, Deshi Shi, Zutao Zhou, Fangping Zhang, Yuncai Xiao.

This repository accompanies the single-cell RNA-seq, wet-lab pharmacology and
chemical-profiling analyses of an APRR herbal formulation in a PRRSV (strain WUH3)
pig challenge model. It is the repository referenced in the *Data availability*
statement of the manuscript.

## Data sources

| Resource | Accession / location | Note |
|---|---|---|
| Raw scRNA-seq reads | Genome Sequence Archive **PRJCA065583** | Controlled access; apply through the Data Access Committee on the accession page |
| Processed Seurat objects (`.rds`) | Not distributed | 11 GB of intermediate objects; see `docs/REPRODUCE.md` to rebuild |
| Derived tables | `tables/` | CSV/TSV, published here |
| Main and supplementary figures | Submitted with the manuscript | Not duplicated here |

**Experimental design.** Three groups (Control, PRRSV, APRR), one animal per
group with cells pooled within group, 24,896 cells after QC (8,488 / 7,794 / 8,614).
Because each group derives from a single animal, group-level comparisons are
descriptive point estimates; statistical inference is performed at the cell level
with false-discovery-rate correction. This limitation is stated in the manuscript.

## Repository layout

```
scripts/                  Final-caliber analysis and figure scripts (R 4.4.2 / Python 3)
  pipeline_r.R            Master pipeline (contains a LEGACY ssGSEA block - see warning below)
  _ssgsea_reproduce_and_fix.R   FINAL ssGSEA: legacy audit + GSVA 2.0.7 curated gene sets
  ssgsea_pathway_activity.R     MSigDB gene-set sensitivity analysis (not adopted)
  nichenet_analysis.R           Literature-curated ligand-target priority scoring
  make_figS9_v2.R               Supplementary Fig. S9 heat map
  _fig4E_func_v2.R              Fig. 4E three-facet functional-program panel
  fig5ij_redraw.py              Wet-lab quantification panels (IL-1beta, p-p65)
  ...                           QC, DEG, CellChat, composition and panel scripts
tables/                   Derived tables used for the reported results
tables/_audit/            Intermediate ssGSEA versions retained for provenance
tables/wet/               Wet-lab summary tables (viral load, phenotype, cytokines, IF)
docs/                     Reproducibility, method notes, software environment
```

## Important: which ssGSEA result is the published one

An audit of the original pipeline found that its pathway-activity block computed
the **inverse** of an enrichment score (mean expression *rank*, where rank 1 =
highest expression, divided by the number of genes) over a 130-gene panel only.
The value was nonetheless labelled "ssGSEA score". The table was reproduced
exactly from the original code, confirming that the stored numbers were correct
but the algorithm was not.

The published scores are from **`_ssgsea_reproduce_and_fix.R`**: GSVA 2.0.7
`ssgseaParam(..., normalize = TRUE)` over the full 17,336-gene transcriptome,
using the same 12 curated gene sets. `tables/ssGSEA_pathway_activity.csv` is the
corrected table used for Supplementary Fig. S9. The legacy table is kept in
`tables/_audit/` as `ssGSEA_pathway_activity_LEGACY_pipeline-rankmean.csv`.
Do **not** use `pipeline_r.R`'s pathway block to reproduce the figure.

The ligand-priority ranking in `nichenet_analysis.R` is a **literature-curated
heuristic**, not the `nichenetr` package model: it scores ligands by the overlap
between a curated ligand-receptor-target dictionary and measured DEGs, with an
explicit directional prior. The manuscript describes it accordingly and treats it
as hypothesis prioritization rather than an unbiased estimate of ligand activity.

## Reproducing

See [`docs/REPRODUCE.md`](docs/REPRODUCE.md). All absolute paths in the released
scripts were rewritten to repository-relative paths; set the working directory
to the repository root before running them. Raw reads must be obtained from GSA
PRJCA065583 first.

## Licence

Code: MIT (`LICENSE`). Derived tables: CC0 1.0 public domain dedication.
If you reuse these data, please cite the manuscript and the GSA accession.

## Contact

Yuncai Xiao (xyc88@mail.hzau.edu.cn) and Fangping Zhang (zfp0632@163.com),
College of Veterinary Medicine, Huazhong Agricultural University, Wuhan, China.
