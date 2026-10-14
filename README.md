# Acne — Molecular Research & Experimental Therapeutics

**Status: exploratory research; no acne cure, drug, or clinically validated prediction has been established.**

This repository collects the October 14, 2026 acne research investigation, including mechanistic analysis of hormonal, sebaceous, follicular, microbial, genetic and immune factors; the FAR/FAR-2 hypothesis; a hypothetical intervention-space screen; and experimental plans for safety and durable remission.

## Latest research — Round 4: actual human acne transcript counts (October 14, 2026)

The repository now includes an independent descriptive analysis of **all 24 public Xenium cell-feature matrices** in NCBI GEO GSE301280 (345,897 recorded cells, 541 matrix features per sample). We identified a public GEO title-count mismatch versus the paper's stated healthy/nonlesional counts and found that restricting expression to marker-detected cells changes apparent FABP5/CRABP2 patterns markedly.

- [Human-data audit with exact group summaries and caveats](research/2026-10-14/ROUND4_GSE301280_PRIMARY_DATA_AUDIT.md)
- [Reproducible GEO downloader and counts analysis](research/2026-10-14/gse301280_small_matrix_audit.py)
- [Derived public data (sample and group CSVs)](research/2026-10-14/gse301280_results)
- [Successful independent GitHub Actions run](https://github.com/f8bsd9yqgp-maker/Acne/actions/runs/38644879207)

**Limits:** This is not a donor-level significance test, a validated cellular-identity analysis, treatment trial, drug discovery or cure. It requires reconciliation with the annotated h5ad and original donor metadata.

## Latest research — Round 3 (October 14, 2026)

The third research phase sharpened the central **FAR-2** question: could a patient-derived human follicle retain an improved *functional* state after treatment withdrawal, and how could a slow return to disease be distinguished from genuine multistability?

- [Round 3 mathematical and biological research audit](research/2026-10-14/ROUND3_STATE_IDENTIFIABILITY_AND_MODEL_FEASIBILITY.md)
- [Reproducible rate-dependent hysteresis diagnostic](research/2026-10-14/state_switch_diagnostic.py)
- [FAR-2 decisive experiment](research/2026-10-14/FAR2_DECISIVE_EXPERIMENT.md)

**Important:** The hysteresis study is a **synthetic generic dynamical-systems exercise**. It does not represent experimental data about acne, demonstrate a lasting follicular reset, or identify a clinically viable molecule. A bug in the original archive's single-target comparator export has now been fixed in the GitHub code and checked by the automated reproduction workflow.

## Research direction

1. Explain why similar hormonal and microbial exposure can lead to different acne severity.
2. Investigate follicular differentiation, sebum lipid quality and inflammation while preserving the skin barrier.
3. Evaluate whether disrupting bacterial virulence amplification (including *C. acnes* GehA) could help without eradicating the normal microbiome.
4. Test whether temporary interventions can produce **drug-free** remission rather than merely transient suppression.
5. Prioritize independently replicated human tissue experiments before any clinical claims.

## Reproducing the conceptual screen

Python 3 and dependencies in `requirements.txt` are needed:

```bash
python -m pip install -r requirements.txt
cd research/2026-10-14
python screen.py
python sensitivity.py
```

The Python screen enumerates **102,631 hypothetical intervention profiles** under assumed biological coefficients. They are **not chemical compounds**, and its scores, simulated improvements and safety penalties are *not validated estimates* of acne outcomes or actual skin safety.

## Research integrity and safety

- Experimental target hypotheses include GATA6/follicular differentiation, FABP5/CRABP2/retinoid routing, lipid oxidation and *C. acnes* GehA.
- Existing studies motivate these targets but do not prove that FAR-2 or any combination prevents human acne or safely provides lasting remission.
- Selectivity, independent donor replication, barrier measurements, tissue penetration, withdrawal and recurrence testing are prerequisites before human trials.
- **Do not try experimental chemical or genetic interventions on your skin.** This repository is a research workspace, not personalized treatment advice.
- The attached original project ZIP in the ChatGPT conversation also contains the full historical reports, notebook, plots and generated model data. Not every large or binary asset is included in this initial GitHub import.

Last research update: **2026-10-14**.
