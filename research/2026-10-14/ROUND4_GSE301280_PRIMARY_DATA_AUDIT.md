# Round 4 — Direct reanalysis of public human acne Xenium count matrices
**Date:** October 14, 2026  
**Data:** NCBI GEO GSE301280, 24 individual `cell_feature_matrix.h5` files (GSM9079610–GSM9079633).  
**Code:** [gse301280_small_matrix_audit.py](gse301280_small_matrix_audit.py).  
**Automated run:** https://github.com/f8bsd9yqgp-maker/Acne/actions/runs/38644879207  
**Status:** All 24 public matrices retrieved and analyzed successfully. This is an independent **descriptive technical replication**, not independent-donor hypothesis testing, causal evidence, proof of a biological attractor, or validation of a treatment.

## 1. Methods and unit of analysis

Every available original Xenium `cell_feature_matrix.h5` was downloaded from NCBI GEO and parsed with `h5py` and the 10x sparse matrix layout; each contains **541 reference features**. Across 24 sections the matrices contain **345,897 recorded cells** (5,384–42,558 per section). We extracted selected transcripts including GATA6, FABP5, CRABP2, KRT5, KRT79 and related differentiation/lipid markers.

For gene `g` and cell `j` we computed `log1p(10000 * count_gj / max(total_panel_counts_j,1))`, then averaged within each section. Tables summarize **medians of section means**, not 345,897 statistically independent observations. This normalization differs from the full original paper's annotated h5ad and is not corrected for tissue processing or donor structure.

We also computed **exploratory marker-positive sensitivity checks** limited to cells with `KRT5>0` or separately `KRT79>0`. These are **not validated cell-type assignments**: a marker's detectability depends on cell state, section quality and read depth. Conditioning on detectable markers could itself introduce selection bias. We did **not** import the authors' `bulk_labels` cell annotations and did not infer donor IDs.

## 2. Metadata discrepancy

| Group | 24 GEO sample titles | Count stated in Durgin et al. 2026 |
|---|---:|---:|
| Healthy | 4 | 5 |
| Nonlesional | 10 | 9 |
| Comedonal | 7 | 7 |
| Pustular | 3 | 3 |

This 1-sample discrepancy between GEO *title labels* and the publication is an **unresolved annotation-reconciliation question**, not an allegation of flawed primary scientific results. Manuscript methods report **11 donors**, but deposited sample accession titles cannot be presumed to provide donor identities; our analysis therefore has no legitimate donor-level p-values.

## 3. First result: whole-section expression is affected by composition

All cells, section medians in log1p counts per 10,000 panel-feature UMIs:

| Tissue (GEO title) | GATA6 | FABP5 | CRABP2 | KRT5 |
|---|---:|---:|---:|---:|
| Healthy (n=4 sections) | 0.542 | 2.282 | 2.193 | 3.264 |
| Nonlesional (n=10) | 0.490 | 2.363 | 1.850 | 2.355 |
| Comedonal (n=7) | 0.206 | 2.335 | 1.483 | 1.954 |
| Pustular (n=3) | 0.181 | 2.707 | 1.593 | 1.696 |

GATA6 is lower in lesion samples, but whole-section KRT5 also falls. The nearly flat nonlesional-to-comedonal FABP5 comparison in all cells does **not** rule out an expression increase within marker-positive cells.

## 4. Second result: marker-gated comparisons are different

Median of per-section means, **only among KRT5-detected cells**:

| Tissue (GEO title) | GATA6 | FABP5 | CRABP2 |
|---|---:|---:|---:|
| Healthy (4) | 0.787 | 3.106 | 3.065 |
| Nonlesional (10) | 0.760 | 3.126 | 2.636 |
| Comedonal (7) | 0.372 | 4.013 | 2.920 |
| Pustular (3) | 0.293 | 5.271 | 3.987 |

Independently restricting to **KRT79-detected cells** yielded:

| Tissue (GEO title) | GATA6 | FABP5 | CRABP2 |
|---|---:|---:|---:|
| Healthy (4) | 0.951 | 3.775 | 3.619 |
| Nonlesional (10) | 0.866 | 3.625 | 3.072 |
| Comedonal (7) | 0.559 | 4.737 | 3.968 |
| Pustular (3) | 0.324 | 6.001 | 4.829 |

**Important finding:** FABP5 appears largely flat comparing nonlesional against comedonal sections when all cells are pooled, but is higher in both marker-positive gates. CRABP2 similarly changes apparent direction relative to the all-cell aggregate. GATA6 stays lower in comedonal and pustular sections even after these two simple gates.

This is **consistent with composition masking or heterogeneity** and gives a strong motivation for genuine cell-type-specific analysis; it is **not formal demonstration of Simpson's paradox**, nor is it proof that FABP5 directly causes acne. Gating thresholds could cause selection effects. The principal publication independently reported spatially localized FABP5 alterations; we should not claim that association as an unprecedented biological discovery.

## 5. Audit completeness and missingness

- All **24/24** H5 matrices downloaded; workflow succeeded.
- Total **345,897 cells**; 541 feature names per matrix.
- Of 20 requested gene symbols, **14 were present**: GATA6, FABP5, CRABP2, KRT5, KRT79, KRT10, PPARG, AWAT2, RARA, RXRA, MKI67, IL1B, SCD, FASN.
- **6 missing panel gene names**: JUN, JUNB, FOS, FOSB, RARB and RXRB. Missing is **not** zero expression. The 541 matrix features and original paper's custom 100-gene assay description need reconciliation with negative controls and feature filtering.
- No inferential claims about significance, response durability, or a drug target are warranted.

## 6. What this does and does not establish

**Established by our reanalysis:** independent retrieval of raw public count files, reproducible per-section gene summaries, a title-label count discrepancy requiring resolution, and strong sensitivity of FABP5/CRABP2 comparisons to a cell-marker-gating choice.

**Not established:** donor-independent robustness; matched cell-type effect sizes; predictive ability for lesion onset; causal links; bistability; an off-treatment reset; treatment safety; a novel acne cure.

**Next scientific gate:** use the original ~758 MB combined **annotated h5ad** and source paper's supplement to recover `bulk_labels`, true donor and sampling metadata, archival-versus-recent flag and the published lesion categories. Re-evaluate GATA6 and FABP5 by verified upper-follicle and sebocyte cell compartments with donor-aware aggregation. If no trustworthy donor mapping is possible, seek clarification before making donor-level statistical inferences.

## Sources

- Public series and section titles: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280
- Durgin et al., *JCI Insight*, February 9, 2026: https://insight.jci.org/articles/view/198021
- Original segmentation repository: https://github.com/joe-durgin/TranscriptFocusedSegmentation

All produced datasets and aggregated results are derived from already-public source data; they are **not patient-specific diagnosis or advice**.
