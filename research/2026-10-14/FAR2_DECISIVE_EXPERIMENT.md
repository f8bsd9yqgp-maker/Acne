# FAR-2 research extension: Can acne-prone human follicles show reversible state switching?

**Date: 2026-10-14. Status: hypothesis, causal-inference plan and literature audit; no biological experiments performed.**

## 1. Decision
The next defensible scientific target is **identifying (or falsifying) a persistent, reversible acne-relevant follicular state under matched external conditions**. This is distinct from merely observing a low GATA6/high FABP5 expression pattern in lesions or a transient benefit of a retinoid.

Operational working hypothesis:
> In a measurable subset of patient-derived human follicular epithelium, the same matched hormonal, microbial and culture environment can sustain two distinguishable functional states of comedogenesis susceptibility after a transient reversible perturbation and confirmed withdrawal of active treatment.

It does **not** predict that all human acne is bistable, or that expression of GATA6 alone provides a cure.

## 2. New literature cross-check and novelty guardrail
- Oulès et al. (2020, *Nat Commun*, doi:10.1038/s41467-020-18784-z): GATA6 supports normal upper pilosebaceous differentiation; several experiments involve immortalized sebocytes and 3D organoids. Not evidence of off-drug reset.
- Durgin et al. (Feb 2026, *JCI Insight*, doi:10.1172/jci.insight.198021): a targeted spatial panel in healthy, nonlesional, comedonal and pustular skin; elevated FABP5/CRABP2, AP-1 activity and reduced retinoid receptor expression in lesional cells. T-5224 lowers FABP5 and decreases pustules in a *mouse folliculitis* model. This does not identify a FABP5 causal contribution to human comedogenesis, nor state switching.
- Balogh et al. (2026, *J Invest Dermatol*, doi:10.1016/j.jid.2025.10.610): cultured breast/abdominal human epidermal keratinocytes can exhibit **C. acnes**-triggered innate immune memory. Not demonstrated for patient-derived facial follicular acne states.
- Bernard (2012, *Exp Dermatol*, doi:10.1111/j.1600-0625.2012.01457.x): proposed a *bistable model of hair growth/dormancy*. This is separate from and is NOT validation of a two-state model of acne susceptibility.
- Lai et al. (2026, *J Invest Dermatol*, doi:10.1016/j.jid.2026.07.015): GehA positively amplifies *C. acnes*/sebum interactions, and neutralization has mechanistic effects. GehA alone was insufficient to induce acne-like pathology in the animal model. A GehA blocker therefore cannot be presumed to cure acne.

**Prior-art limit:** This is not a claim of first-ever epithelial memory, first-ever GATA6 manipulation, or first-ever mathematical bistability in hair follicles. A focused originality question is whether *acne-relevant upper-follicle histopathology* can be reversibly switched between persistent states with measured common inputs. Even that claim requires broader publication and patent searches before assertion.

## 3. A mathematical obstacle: apparent persistence cannot prove bistability
Consider two dynamical accounts of a transiently normalized follicular readout `x(t)` measured after treatment has stopped.

- **Unique equilibrium with slow relaxation:** `dx/dt = -(x-x_disease)/tau`, so `x(t)=x_disease + (x(0)-x_disease)*exp(-t/tau)`. With sufficiently large tau, this appears stable over any finite observation period.
- **Candidate bistability:** `dx/dt = x - x^3 + h - u(t)` (toy model only). At constant h with `|h| < 2/(3*sqrt(3))`, there may be two stable state equilibria separated by an unstable equilibrium.

**Identifiability statement:** Finite-duration persistence alone cannot distinguish an arbitrarily slowly relaxing monostable model from a truly bistable one at arbitrary precision. For any finite time window `T` and tolerance `epsilon>0`, choose `tau` large enough that `1-exp(-T/tau)<epsilon`. The readout scarcely changes, although the only long-time equilibrium is the acne-prone state.

Consequently, a single post-washout endpoint—even with cleared drug—is insufficient to establish an alternative attractor.

## 4. More discriminating approach: common-input switching challenge
Design within-donor, *not merely between-donor*, matched follicular-model experiments:

**Stage A: Validity and baseline**
1. Independently replicate abnormal follicle differentiation in a human acne-derived experimental system; require morphologic microcomedone or validated infundibular occlusion readouts, not just biomarker ratios.
2. Assess whether untreated sister cultures remain sufficiently stable under same medium, hormonal and microbial stimuli.
3. Measure GATA6, FABP5 and CRABP2 protein, RAR/RXRA competence, RAR functional outputs, differentiation markers and lipids. Comedonal and pustular forms must be considered separately.
4. Check donor, anatomic location, recent treatment, donor sex/hormonal context, passage, culture-media drift, cell-type mixtures and cell viability.

**Stage B: Mechanistic target isolation**
Randomized within-donor arms: vehicle, reference retinoid, orthogonal FABP5-specific perturbations, matched specificity rescue, CRABP2 perturbation, and appropriately designed GATA6 perturbation/rescue. Separately compare non-specific AP-1 inhibition. Do not attribute AP-1 benefit to FABP5 without genetic rescue and pathway readouts.
- Readout must be **morphology/occlusion** plus molecular biomarkers; reduction of a transcript is not sufficient.
- Knockdown and rescue delivery controls must exclude transient toxicity, off-target effects and paracrine confounding.
- Changes to GATA6 must be evaluated for TGF-beta dysregulation, normal barrier differentiation and immune function.

**Stage C: Persistence versus hidden pharmacology**
Randomize *matching models from the same donors* into pulse and sham conditions. Stop the intervention and verify the agent is cleared using appropriate analytical and functional assays. Assess molecular pathway reversion, cell turnover and maturation.
- Same standardized environment is maintained in both arms after treatment.
- Compare repeated, longitudinal morphology, not merely a single endpoint.
- Culture duration must be biologically sufficient. Existing pluripotent-stem-cell skin organoids can be maintained up to ~150 days in some protocols, but resemble fetal tissue and have not been validated as faithful models of acne comedogenesis (doi:10.1038/s41596-022-00681-y).
- If a model cannot remain stable long enough, classify the result as indeterminate rather than a reset.

**Stage D: Reversible switch and hysteresis**
At identical environmental input settings, compare outcomes after two distinct preparation histories (pulse versus untreated). Test state separation across increasingly long dwell times; reverse the perturbation where biologically justified and ethically approved.
- Evidence favoring two-state behavior: sustained morphological separation under the *same* inputs; within-donor reproducible switching in both directions with predictable pulse thresholds, convergence within each state and preservation of cell function.
- Control for delayed signaling by varying pulse duration, recovery time and environmental ramp speed. Rate-dependent lag is **not** evidence of equilibrium hysteresis.
- Reject claims of persistent biological attractors if states converge after enough time, if residual target engagement explains differences, if differentiation merely reflects cell-population selection, or if model viability declines.
- Even positive model findings support *model-specific multistability*, not a permanent human cure.

## 5. Experimental decision tree
- **If perturbation improves biomarkers but not follicular architecture:** target mechanism is insufficient; do not proceed with cure framing.
- **If a benefit disappears at effective agent clearance:** conventional symptomatic effect, not FAR-2 drug-free reset.
- **If sustained benefits remain but only one switching direction works:** promising durability, but not proof of bistability; pursue alternative mechanisms including delayed maturation and population selection.
- **If reversible switching persists across independent donors with common inputs and no barrier/immune cost:** justify expert collaboration, replication and long-term follow-up.
- **If only a donor subgroup responds:** prioritize predictors and subgroup-specific hypotheses; do not claim universal mechanism.

## 6. Study design and statistics requirements
- Independent unit is **the donor**, with condition randomization nested within donor; technical cells/follicles are not independent biological replicates.
- Investigators preregister a primary morphological endpoint, time windows, rescue effects, meaningful margin, tolerable barrier harm and dropout criteria.
- Collect longitudinal raw readouts and agent-clearance data, masked outcome scoring and donor-level mixed-effects analyses.
- Test both spontaneous recovery and return under endocrine/microbial challenge. Be aware artificial challenges may not reproduce patient acne.
- Do not promise a sample size from the model. Use an independent pilot to estimate donor heterogeneity and variability, then pre-specify a power or precision design.
- A failed model is uninterpretable when the negative control also loses follicular architecture or model viability.

## 7. Priority and scope
**Immediate:** Analyze published public datasets for donor-level and lesion-subtype-specific upstream signals, alongside feasibility of a human follicle model. This does not require inventing a new drug.
**Next:** Ask an experienced human pilosebaceous-unit/acne researcher to critique Stage A/B feasibility, morphology measures, time scale, and donor-level sample-size requirements.
**Later:** Only after target causality and robust in-vitro recovery: target-specific delivery engineering and safety/clinical evaluation.

No new wet-lab data, donor tissue, drug molecules or validated target-directed interventions were produced by this report. Existing acne treatment and human research ethics remain the appropriate standards.

## Sources verified / links
1. https://www.nature.com/articles/s41467-020-18784-z
2. https://insight.jci.org/articles/view/198021
3. https://pubmed.ncbi.nlm.nih.gov/41242652/
4. https://doi.org/10.1111/j.1600-0625.2012.01457.x
5. https://pubmed.ncbi.nlm.nih.gov/42526582/
6. https://pubmed.ncbi.nlm.nih.gov/35322210/
