# Software environment

Measured on the analysis workstation (Windows, R 4.4.2):

| Package | Version |
|---|---|
| R | 4.4.2 |
| Seurat | 5.5.1 |
| SeuratObject | 5.4.0 |
| GSVA | 2.0.7 |
| ggplot2 | 4.0.3 |
| dplyr | 1.2.1 |
| Matrix | 1.7.5 |
| patchwork | 1.3.2 |
| scDblFinder | 1.20.2 |
| SingleCellExperiment | 1.28.1 |
| nichenetr | 2.2.1.1 (installed; **not used** for the reported ranking) |

CellChat version, library kits, sequencing platform, antibody/catalogue numbers
and all animal-procedure details are reported in the manuscript Methods section;
refer to that section as the authoritative description of the analysis as
performed. Python scripts require Python 3 with pandas, matplotlib and Pillow.

Sequencing: Illumina NovaSeq 6000, paired-end; 10x Genomics Chromium Next GEM
Single Cell 3' Reagent Kits v3.1, target recovery 8,000 cells per sample.
Reference genome: Sscrofa11.1; alignment/quantification with Cell Ranger v7.0.1.
