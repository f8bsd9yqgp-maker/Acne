# Acne — Molecular Research & Experimental Therapeutics

**Status: exploratory research; no acne cure, drug, or clinically validated prediction has been established.**

This repository collects the October 14, 2026 acne research investigation, including mechanistic analysis of hormonal, sebaceous, follicular, microbial, genetic and immune factors; the FAR/FAR-2 hypothesis; a hypothetical intervention-space screen; and experimental plans for safety and durable remission.

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
