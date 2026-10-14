# Acne conceptual screening model — 2026-10-14

This is **hypothesis-generation research only**. The controls represent possible therapeutic *mechanisms*, not molecules or clinically tested treatments.

- `screen.py` performs an exhaustive deterministic grid of 102,631 conceptual target interventions under 4 synthetic phenotypes and 4 hypothetical biological coefficient sets. It generates all CSV outputs, two PNG charts, and a plain-text summary.
- `sensitivity.py` reranks the model's results with alternate hand-picked scoring objectives and writes `objective_sensitivity.csv`.
- Run `python -m pip install -r requirements.txt` in the repository root; then `cd research/2026-10-14 && python screen.py && python sensitivity.py`.

**Interpretation:** The percentages are unvalidated reductions in a mathematical risk proxy—not measured lesion-count reductions, clinical efficacy, docking scores, ADMET predictions or evidence of skin safety. The 102,631 profiles are not candidate drugs. All modeled coefficients and guardrails are assumed rather than fitted to clinical outcomes.

The full historical reports, original notebook, images and precomputed CSVs are included in the separately supplied research ZIP.
