# Method notes and audit trail

## ssGSEA: why there are several tables

| File | Origin | Used in manuscript |
|---|---|---|
| `tables/ssGSEA_pathway_activity.csv` | `_ssgsea_reproduce_and_fix.R`, GSVA 2.0.7, 12 curated gene sets, full transcriptome | **YES** (Fig. S9) |
| `tables/_audit/ssGSEA_pathway_activity_LEGACY_pipeline-rankmean.csv` | `pipeline_r.R` rank-mean over 130 genes | No |
| `tables/_audit/ssGSEA_pathway_activity_reproduced.csv` | Exact byte-level reproduction of the legacy table | Audit only |
| `tables/_audit/ssGSEA_pathway_activity_rebuilt.csv` | GSVA with MSigDB Hallmark/Reactome/GO sets | Sensitivity analysis only |
| `tables/_audit/ssGSEA_pathway_activity_gsva_curated.csv` | GSVA with the 12 curated sets (pre-normalisation working copy) | Audit only |

The legacy algorithm ranks genes within each cell (rank 1 = highest expression),
then averages the ranks of gene-set members and divides by the number of genes.
A *lower* rank therefore means *higher* expression, so a high "score" indicated
low pathway expression - the opposite of the stated interpretation. Its gene
universe was also restricted to 130 marker / ligand-receptor / pathway /
checkpoint genes rather than the full transcriptome.

Reproduction check: the legacy code re-run on the current 3-group object
reproduced the archived table with a maximum absolute difference of 0.00000000
and Spearman rho = 1.0000 for all 12 pathways. The stored numbers were correct;
the computation was not a valid ssGSEA score.

With the corrected GSVA scores the direction of the interferon rows agrees with
the manuscript (PRRSV suppresses IFN-alpha / antiviral programmes; APRR partly
restores them), whereas the legacy table implied the opposite. No main-text
conclusion depended on the ssGSEA values; only Supplementary Fig. S9 was
affected and was regenerated.

## Ligand-target priority score

`nichenet_analysis.R` does not run the `nichenetr` package. It computes

```
priority = log1p(sender expression)
         x (n_up - 0.3*n_down + 0.5*n_recovery + 0.5*max(R, 0))
         x d
```

where `n_up`, `n_down` and `n_recovery` count overlaps between a
literature-curated ligand-receptor-target dictionary and measured DEG sets,
`R` is the receptor-activity term and `d` is a hand-assigned directional prior
(e.g. IL1B = +1, IL10 = -1). Because `d` is a prior, the top rank of IL1B is
partly specified by the model. The manuscript therefore calls this a
literature-curated ligand-target scoring / hypothesis prioritization and does
not treat it as independent mechanistic evidence.

## T-cell functional panel

The original two-facet classification (`Cytotoxicity/effector` vs
`Exhaustion/checkpoint`) was a binary fallback: every gene not in a list of six
cytotoxicity genes was placed in the "exhaustion" facet, which incorrectly
included CCR7, TCF7, CD44, IL7R and CD69. `_fig4E_func_v2.R` uses three
facets (`Cytotoxicity/effector`, `Checkpoint/inhibitory`,
`Progenitor/activation`) with an explicit gene-to-facet map and a `stopifnot`
guard. The manuscript text reports per-gene directions rather than facet labels.

## Statistics

Two-group comparisons of viral load and clinical endpoints use unpaired
Student's t-tests, except serum IL-1beta and lung p-p65 (Welch's t-test).
Group-level single-cell comparisons are descriptive; inference is at cell level
with FDR correction. `_audit_stats_recompute.R` independently recomputes the
in-vivo statistics from the raw wet-lab tables.
