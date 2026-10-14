# Round 5 — Full annotated human spatial transcriptomics reanalysis
**Completed 2026-10-14.** Public data: **GSE301280**; original research Durgin et al. (2026).  
**Status:** Independently retrieved the deposited **757.9 MB h5ad**; inspected original cell-type/condition annotations; calculated expression by tissue section and annotated cell type. No independent donor-level significance, perturbation, longitudinal recurrence, drug safety or acne cure is established.

## 1. Provenance and pipeline

- Public annotated matrix: `GSE301280_anndata_global_6_29_25.h5ad` at NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280).
- Code: `gse301280_annotated_h5ad_audit.py`, `annotated_sensitivity_report.py`. GitHub Actions complete: https://github.com/f8bsd9yqgp-maker/Acne/actions/runs/38645873697
- `ad.read_h5ad(backed="r")`, separate 3,000-cell block reads, 14 preregistered-style **exploratory target genes** (not pre-registered before observation), and grouping on deposited `bulk_labels`, `prefix_healthy_label`, `sample` and `new` metadata.
- **258,965 retained cells**, **100 genes**, **15** categorical `bulk_labels` cell types. Public unannotated feature matrices examined in Round 4 contained 345,897 cells and 541 total feature names per sample: the data representations and filtering are **not interchangeable**.
- Gene expression is read directly from the deposited annotated `X`, already normalized/log-transformed according to GEO's processing description. All values below are **median across tissue-section means** of X, *not* arbitrary transformed fold-change, cell-level inferential tests or clinical effect sizes.

## 2. The sample-classification discrepancy

The annotation `prefix_healthy_label` yields the **same** count breakdown the article describes:

| Annotated state | sections |
|---|---:|
| Healthy | 5 |
| Acne Non-lesional | 9 |
| Acne Comedo | 7 |
| Acne Pustule | 3 |

Public GEO *sample accession titles* instead appeared to yield 4 healthy / 10 nonlesional / 7 comedonal / 3 pustular, as Round 4 reported. The mismatch is **reconciled at aggregate count level** by using the authors' disease labels; exactly which accession title maps to any reclassified internally numbered sample remains to be verified. Do not assume paper mislabeled a sample.

## 3. Collection-era batch structure and missing donor IDs

`new` splits **18 recent and 6 archival sections**, consistent with the publication's archival description. The six archival sections are **3 nonlesional + 3 comedonal**; all five healthy and three pustular sections are recent. This means lesion category and collection era are partially confounded. Analyses restricted to recent sections provide an important, but not sufficient, sensitivity check.

The annotated `obs` metadata has `sample` (24 unique), `specimen` (19 unique), `specimen_id` (24 unique), and other technical fields, but **no verified patient/donor identifier**. The paper reports **11 donors**. None of these numbered technical fields may responsibly be treated as 11 independent donors without documented mapping. Thus we do not report p-values, false-discovery rates, cross-donor effects, or clinical prediction.

## 4. Key descriptive results in annotated cells

Median tissue-section means from h5ad X, disease label sourced from `prefix_healthy_label`.

| Cell population | Gene | Healthy (5) | Nonlesional (9) | Comedo (7) | Pustule (3) |
|---|---|---:|---:|---:|---:|
| Hair-follicle KC | GATA6 | 0.508 | 0.295 | 0.230 | 0.033 |
| Hair-follicle KC | FABP5 | 0.929 | 1.075 | 1.157 | 1.653 |
| Basal KC | GATA6 | 0.059 | 0.053 | 0.010 | 0.008 |
| Basal KC | FABP5 | 0.291 | 0.298 | 0.444 | 0.662 |
| Transitional KC | GATA6 | 0.118 | 0.111 | 0.045 | 0.003 |
| Transitional KC | FABP5 | 0.690 | 0.736 | 1.109 | 1.566 |

This is a **consistent descriptive pattern** in multiple keratinocyte labels: lesion sections tend to have lower GATA6 and higher FABP5, matching the direction of the published study's emphasis. Because these are cross-sectional, observational, normalized count data, the pattern **neither identifies causality nor establishes an independent discovery**.

### Sensitivity restricted to recently collected sections

| Cell population | Gene | Healthy (5) | Nonlesional (6) | Comedo (4) | Pustule (3) |
|---|---|---:|---:|---:|---:|
| Hair-follicle KC | GATA6 | 0.508 | 0.429 | 0.119 | 0.033 |
| Hair-follicle KC | FABP5 | 0.929 | 1.216 | 1.624 | 1.653 |
| Basal KC | GATA6 | 0.059 | 0.046 | 0.006 | 0.008 |
| Basal KC | FABP5 | 0.291 | 0.300 | 0.394 | 0.662 |
| Transitional KC | GATA6 | 0.118 | 0.075 | 0.003 | 0.003 |
| Transitional KC | FABP5 | 0.690 | 0.655 | 1.413 | 1.566 |

Within recent samples, basal and transitional keratinocyte comparisons retain a similar directional pattern. Sample counts are **small**, and donor dependence remains unknown.

## 5. Sensitivity to tiny annotated cell populations

Do not interpret every cell-type comparison equally. For HF KC, per-section annotated counts were:

- Healthy: 82–609; nonlesional: 71–983; comedo: **7–425** (two comedo sections had only 7 or 17 HF KCs); pustule: 419–1,618.
- With a **minimum of 100 cells** per sample: HF KC coverage becomes 4 healthy, 7 nonlesional, 5 comedo, and 3 pustule sections. The all-era comedonal GATA6/FABP5 medians remain about **0.230/1.157**; however, *recent comedo* HF KCs have only **2** qualifying sections, so their trends are fragile.
- Basal KC counts are at least **548** cells in every analyzed section, and Transitional KC at least **214**. These are stronger **coverage-supported** descriptive compartments.
- **Sebocytes are too sparse in several lesions for robust section comparisons.** Comedo section counts: `1,1,3,40,449,784,954`; pustule: `4,12,339`. Only 3/7 comedo and 1/3 pustule sections have ≥100 sebocytes. Claims of a smooth sebocyte trend would be misleading without larger, better-resolved samples.

All full tables and exclusion rules can be reproduced with `annotated_sensitivity_report.py` from the stored section-by-cell-type CSV.

## 6. Interpretation relevant to FAR-2

- **What strengthened:** biologically interpretable, original annotated data support revisiting **follicular/keratinocyte retinoid transport and differentiation** and suggest the FABP5/GATA6 associations are not *solely* an artifact of pooling entirely unrelated cell populations. This is *supportive*, not a mechanism proof.
- **What did not strengthen:** the dataset is cross-sectional (no baseline-to-lesion follow-up, no perturbation, no washout), so it provides **zero direct evidence** for FAR-2's distinct long-lived attractors or any permanent acne reset.
- **What weakens naive assertions:** poor sebocyte representation in many active lesions and category/archival differences. Low GATA6 expression in some subtypes approaches a detection floor and should be checked with independent protein/spatial assays.
- **What is not new:** the original researchers already reported elevated FABP5 and differentiation/retinoid signaling abnormalities. Our original work here is **an independently executable computational replication plus data-quality and cell-coverage audit**, not first discovery of those biology signals.

## 7. Highest-value next biological step

1. Obtain the actual donor/section crosswalk (possibly from original study investigators or published supplement) and verify whether replicates are within the same donor; apply donor-level aggregation/paired contrasts when justified.
2. Test whether altered FABP5 and GATA6 are concentrated in distinct upper-follicular compartments on spatial coordinates, with independent validation for sparse GATA6 expression; control for sample processing era.
3. Require mechanistic perturbation + genetic rescue in reliable acne-relevant human follicular models to establish whether changing the pathway alters **follicular plugging** while maintaining barrier integrity.
4. Only afterward evaluate post-washout durability and state-switching; cross-sectional expression patterns cannot identify bistability.

## Reproducibility links
- https://github.com/f8bsd9yqgp-maker/Acne/blob/main/research/2026-10-14/gse301280_annotated_h5ad_audit.py
- https://github.com/f8bsd9yqgp-maker/Acne/tree/main/research/2026-10-14/gse301280_annotated_results
- https://github.com/f8bsd9yqgp-maker/Acne/blob/main/research/2026-10-14/annotated_sensitivity_report.py
- Primary data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280
- Original publication: https://insight.jci.org/articles/view/198021
