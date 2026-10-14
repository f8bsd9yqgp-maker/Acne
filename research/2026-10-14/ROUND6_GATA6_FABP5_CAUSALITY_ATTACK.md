# Round 6 — Are GATA6↓ and FABP5↑ drivers, consequences, or feedback mediators in acne?

**Date:** 2026-10-14  
**Scope:** Critical original-literature audit, independently recomputed human-data sensitivity, formal causal-identification obstruction, and falsifiable proposed experiment.  
**Status:** No causal genetic perturbation of the public acne cohort and no new wet-lab findings. No new acne treatment or cure established.

## Executive decision

**Prioritize testing GATA6's necessity/sufficiency for morphologic follicular normalization, while directly isolating FABP5's role as a potential modulator (or compensatory consequence), not assuming FABP5 suppression must help.** Strong mechanistic evidence already supports a GATA6 role in human follicular differentiation and retinoid response, but full human acne-prevention and durability have not been shown. FABP5 has experimental evidence for retinoid signaling in non-acne cell systems, but acne-specific causality and direction of benefit remain unresolved.

## 1. Perturbation evidence ledger

| Target | Evidence experimentally established | Not established |
|---|---|---|
| **GATA6** | Mouse epithelial Gata6 loss dilates infundibulum/ducts and causes keratotic changes (doi:10.1111/exd.13757). Human follicle explant siRNA knockdown changes infundibular IVL staining, while GATA6 overexpression changes keratinocyte and sebocyte differentiation; knockdown partially reverses retinoic-acid effects in immortalized sebocytes (Oulès et al. 2020, doi:10.1038/s41467-020-18784-z). | Successful reversal/prevention of *human acne comedones*, off-treatment durable remission, patient-level specificity and safety of targeted GATA6 modulation. Mouse follicles are not a faithful human acne model. |
| **FABP5** | HaCaT studies in 2007 showed that FABP5 RNA interference alters retinoid/PPAR transcriptional output. 2015 experiments reported that fatty-acid ligands alter retinoid receptor activation dependent on FABP5; biology depends on ligand context. Durgin et al. 2026: calcium induced FABP5 in human N/TERT keratinocytes; AP-1 inhibitor T-5224 lowered FABP5 mRNA/protein; same inhibitor reduced murine PMA/high-fat-diet neutrophilic folliculitis pustules by **44%**. | Specific FABP5 intervention that prevents/reverses human comedogenesis, mediation of the in-vivo AP-1 inhibitor benefit through FABP5, drug-free follicular state resetting, or a guarantee that reducing FABP5 improves retinoid signaling in all relevant tissues. |

**Important counterevidence / model dependence:** A 2011 HaCaT keratinocyte study found that adding retinoic acid increased an RAR response but did not increase measured PPARβ/δ targets even at high FABP5:CRABP-II, arguing against a universal, receptor-ratio-only rule (PMID:21843636). A 2015 study showed fatty-acid context shifts RA responses and FABP5 depletion can change both RAR and PPAR signaling (doi:10.1038/ncomms9794). Hence 'FABP5 always diverts retinoic acid away from RAR' is too strong to use as a therapy assumption. Measure RAR *and* PPAR-β/δ functional outputs at relevant lipid conditions and multiple target levels.

**Known literature novelty collision:** The 2020 GATA6 work specifically proposed a follicular fate switch; a 2026 bioRxiv preprint (doi:10.64898/2026.07.09.737640) describes a putative microcomedone-like cell population in nonlesional acne skin and a data-driven comedone-switch framework. Our project must not claim to have invented GATA6 control, FABP5 retinoid effects, or the comedone-switch concept.

## 2. Recomputed human-data sensitivity: GSE301280

Input: independently processed original annotated 758 MB h5ad from 24 tissue sections, using deposited `bulk_labels` and `prefix_healthy_label`. This round adds `round6_causal_signature_robustness.py` and `round6_section_level_robustness.csv` with all group contrasts, collection-era exclusions, ≥100-cell per section filters, leave-one-section-out median deltas and **descriptive** Cliff's delta.

### Comedonal vs nonlesional, recent sections only, each cell type ≥100 cells

| Cell type | Gene | Comedonal sections | Nonlesional sections | Median difference in annotated X | Cliff's delta | Leave-one-section-out difference range |
|---|---|---:|---:|---:|---:|---|
| Basal KC | GATA6 | 4 | 6 | -0.0404 | -0.667 | -0.0450 to -0.0359 |
| Basal KC | FABP5 | 4 | 6 | +0.0930 | +0.667 | +0.0503 to +0.1358 |
| Transitional KC | GATA6 | 4 | 6 | -0.0723 | -0.583 | -0.1105 to -0.0342 |
| Transitional KC | FABP5 | 4 | 6 | +0.7581 | +1.000 | +0.6768 to +0.8394 |
| Hair-follicle KC | GATA6 | 2 | 4 | -0.2568 | -1.000 | -0.4338 to -0.0798 |
| Hair-follicle KC | FABP5 | 2 | 4 | +0.3587 | +1.000 | +0.1725 to +0.5448 |

Each direction remains fixed after deletion of *one section*, but this **does not account for linked donors**, age mismatch, selection or confounding. HF KC comparisons are especially weak: only 2 comedonal sections meet the ≥100-cell recent criterion; do not treat these effect magnitudes as robust evidence across people.

A Cliff's delta of +1 in a sample-section contrast means every **section mean** in one group exceeds those in the other group in this particular subset. It does not imply perfect individual-level discrimination, a biomarker accuracy estimate or causation. Uncertainty intervals and p-values cannot be interpreted without true donor IDs and a design-level analysis.

### Why nonlesional samples cannot prove 'before lesion'
They were sampled cross-sectionally from patients with acne, and 'nonlesional' denotes the state of a sampled area at biopsy, **not** documented future lesion emergence. Subclinical microcomedones, inflammation, systemic predisposition, age difference and differences in cell identity may already exist. The four lesion categories are not an experimentally observed sequence of the same follicle over time.

In our sample summaries, healthy versus nonlesional basal/transitional keratinocyte GATA6/FABP5 differences are generally smaller than nonlesional versus lesion differences, compatible with **both** an early trigger amplified by disease and a reaction to established pathology. They are not enough to choose direction.

## 3. Formal obstruction: cross-sectional gene data do not identify causal orientation

Let G, F, and L represent GATA6, FABP5 and acne-lesion status. If every observed joint probability is positive, the same joint distribution can be written as:

`P(G,F,L)=P(G)P(F|G)P(L|G,F)`

or

`P(G,F,L)=P(L)P(F|L)P(G|F,L)`.

Both factorizations correspond to acyclic causal models with the *same observational likelihood* but different arrows. Unobserved hormones, differentiation, microbiota, inflammation and donor effects create additional non-identifiability. More transcripts, even measured perfectly, cannot select these graphs on their own. Experimentally timed interventions, functional endpoints and relevant assumptions are required.

Plausible competing mechanisms:

1. **Primary ductal-identity driver:** diminished GATA6 impairs ductal differentiation and *then* promotes comedogenesis; inflammatory/differentiation signaling increases FABP5 as a secondary response.
2. **Retinoid-routing driver:** altered lipid/AP-1/FABP5 changes receptor activity and differentiation, lowering GATA6 later.
3. **Common-cause/reactive:** hormones/sebum/stress cause plugging first; inflammation and cell-state shifts both suppress GATA6 and increase FABP5.
4. **Dual driver or positive feedback:** GATA6 and FABP5 influence distinct processes, with established plugging further exacerbating both (state persistence possible but unproven).

None is identified by the public h5ad.

## 4. Discriminating experimental design for specialist human tissue labs

**Eligibility first:** A validated organotypic *human upper pilosebaceous* system with reliable de-novo infundibular microcomedone/occlusion morphology, physiological sebum/keratinization and reproducible viability is necessary. Short-lived explants that show only expression changes cannot adjudicate prevention or stable remission.

**A. Separate sufficiency/necessity and direction**

Within-donor experimental units randomized to:
- Negative and non-targeting controls, matched intervention control, and established appropriate acne-relevant retinoid comparator.
- GATA6 reduction and restoration toward a physiologic follicular range, each with independent orthogonal approaches and specificity rescue.
- FABP5 reduction and restoration toward physiologic range, again with orthogonal approaches and resistant/add-back rescue.
- A preregistered **2×2 factorial** of controlled GATA6 and FABP5 states, plus independent manipulations of common upstream lipid/stress/AP-1 drivers.
- Intervention **before** an approved comedogenic challenge versus **after** a measurable microcomedone phenotype has developed, under otherwise matched conditions.

Blind morphology scoring. Do not equate a change in KRT79/KRT10/FABP5 transcripts with a cleared follicle. Primary endpoint: objectively measured infundibular plugging/microcomedone incidence or burden in a model independently shown to reproduce that feature. Co-endpoints: barrier integrity, tissue damage, cell proliferation/lineage identity, secretion/lipid quality, IL-1/chemokine physiology, GATA6 protein, FABP5 protein, functional RAR/RXRA and PPARβ/δ readouts, cell population composition and viability.

**B. Test mediation rather than correlation**

- If AP-1 suppression improves morphology, but direct FABP5 suppression does not, AP-1 benefits likely arise through other targets; FABP5 was potentially a passenger.
- If directly manipulating FABP5 changes validated functional RAR output *and* occlusion, and restoring FABP5 selectively reverses both, that supports a FABP5-mediated causal effect (conditional on specificity).
- If adjusting GATA6 alters occlusion even when FABP5 is held near its baseline range, this supports a partially distinct GATA6 pathway.
- If raising/restoring GATA6 improves morphology but FABP5 forcing reverses that benefit, investigate epistasis and the underlying receptor pathways. Reverse contrasts test feedback.
- If perturbations alter only stress markers, drug metabolism, or immune-cell infiltration but not ductal plugging, distinguish downstream anti-inflammatory benefit from prevention of comedogenesis.
- If the effect occurs only after disease induction and not with prophylactic perturbation, this favors downstream response or disease-stage-specific action (subject to model adequacy).

**C. Factorial interaction contrast**

With outcome `Y` measuring adverse occlusion and intervention labels `G,F ∈ {0,1}`, estimate
`Interaction = E[Y|do(G=1,F=1)] - E[Y|do(G=1,F=0)] - E[Y|do(G=0,F=1)] + E[Y|do(G=0,F=0)]`,
all at the same experimental history, matching input and time. A nonzero estimate implies interaction at the chosen endpoint scale, **not** automatically direct gene–gene regulation or durable cure. Report donor-level intervals and heterogeneity; no sample-size claim without pilot data.

**D. Only if functional effects replicate: test off-treatment durability**

Confirm active intervention and drug pharmacology have cleared. Maintain identical input conditions and check for repeated, viable phenotypic recovery versus slow relaxation; use reversible challenges and ramp-rate checks from Rounds 2–3. No lasting remission or bistability claim can be made on a short washout or a molecular-only readout.

## 5. What would falsify the therapeutic hypothesis?

- No morphological benefit despite verified, on-target GATA6 or FABP5 perturbation in viable valid follicle models.
- Rescue/add-back fails to reproduce reversal of phenotype, suggesting off-target or model confounding.
- Benefits require damaging the barrier, innate immune defenses or normal sebum function.
- FABP5 reduction redirects receptor signaling in an adverse direction in relevant fatty-acid environments.
- Supposed durable benefit disappears after clearance or one tissue turnover; 'reset' claim rejected.
- Directional contrasts vanish in independent donors or when matched for age, microenvironment and lesion history.

## 6. Next actionable work

1. Recover **24 section ↔ 11 donor** mapping from the supplemental clinical table or study authors. Until then, no donor-level p-values.
2. Check spatial localization in intact ductal compartments; the annotated 100-gene panel has limited pathway readouts (e.g. JUN/FOS absent), so do not infer AP-1 flux from our transcript counts alone.
3. Request expert critique of which **human upper-follicle model** can reproduce a new microcomedone without losing viability.
4. Conduct pre-registered, independent-donor lab perturbation/rescue experiments. Do not jump directly to self-experimentation, compound use or human clinical intervention.

## Primary source links
- Oulès et al. (2020) human explant GATA6 knockdown + retinoid mediation: https://www.nature.com/articles/s41467-020-18784-z
- Gata6 knockout mouse upper duct anatomy (2018): https://pubmed.ncbi.nlm.nih.gov/30033638/
- Durgin et al. (2026) AP-1, FABP5, human N/TERT and mouse folliculitis: https://content.jci.org/articles/view/198021/
- Early HaCaT retinoid receptor FABP5 work (2007): https://pmc.ncbi.nlm.nih.gov/articles/PMC1948722/
- HaCaT counterexample to universal ratio rule (2011): https://pubmed.ncbi.nlm.nih.gov/21843636/
- Fatty-acid ligand dependence (2015): https://www.nature.com/articles/ncomms9794
- 2026 preprint comedone switch: https://www.biorxiv.org/content/10.64898/2026.07.09.737640v1
- Raw GSE301280 and annotation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280

**Responsible interpretation:** The meaningful advance here is prioritization and design quality: existing GATA6 causality evidence, FABP5 mechanistic ambiguity, a robust descriptive cellular association, and the precise *experimental manipulations that would distinguish competing theories*. None constitutes validation of a new acne cure.
