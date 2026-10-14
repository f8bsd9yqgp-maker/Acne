# Round 7 — Recovered donor structure, exact 720-assignment bounds, and a feasible human infundibular assay

**Date:** 2026-10-14  
**Research status:** Methodological / literature advance; no new human experiments, chemical candidates, proof of acne remission or clinical cure.

## 1. The missing source recovered

The primary 2026 paper's **Supplemental Table S1** is freely available (Supplemental data PDF, page 8 / 1-indexed):  
https://insight.jci.org/articles/view/198021/sd/pdf/render/1

It explicitly describes samples from 11 *Xenium* donors, with *repeat samples* within patients. These are not independent biological replicates.

| Source table anonymous donor | Recent/archival | Nonlesional (N) | Comedonal (C) | Pustular (P) | Healthy (H) | Total Xenium sections |
|---|---|---:|---:|---:|---:|---:|
| Patient 7 | archival | 1 | 3 | 0 | 0 | 4 |
| Patient 8 | archival | 1 | 0 | 0 | 0 | 1 |
| Patient 9 | archival | 1 | 0 | 0 | 0 | 1 |
| Patient 10 | recent | 2 | 0 | 1 | 0 | 3 |
| Patient 11 | recent | 1 | 0 | 1 | 0 | 2 |
| Patient 12 | recent | 2 | 1 | 1 | 0 | 4 |
| Patient 13 | recent | 1 | 3 | 0 | 0 | 4 |
| Patient 14 | recent | 0 | 0 | 0 | 1 | 1 |
| Patient 15 | recent | 0 | 0 | 0 | 2 | 2 |
| Patient 16 | recent | 0 | 0 | 0 | 1 | 1 |
| Patient 17 | recent | 0 | 0 | 0 | 1 | 1 |
| **Total** | | **9** | **7** | **3** | **5** | **24** |

This resolves the **donor-level category count structure**, but the annotated h5ad still does not include a verified cell-to-donor or section-to-donor ID. Anonymous patient 12 and 13 labels are from the published supplemental document, not identification of people.

The source table records biopsy sites differing within some patients (e.g., face vs back). Donor matching alone therefore does not eliminate *anatomic site* confounding.

## 2. An original analytical workaround: exact partial identification

The original public h5ad defines 24 samples and 15 cell-type labels. Six *recent* samples are N, four are C, three are P. Patients 12 and 13 are the **only recent donors with both N and C samples**; our question is the mean of two within-person differences:
- Patient 12: `mean(1 C) - mean(2 N)`;
- Patient 13: `mean(3 C) - mean(1 N)`.

We do **not** know which individual N/C section belongs to which of patients 10, 11, 12 and 13. Rather than guess, enumerate all category-compatible allocations:
- Assign one of four C sections to patient 12; patient 13 gets the remaining three: **4** possibilities.
- Assign two of six N to patient 10: **15** possibilities.
- Assign one of the four remaining N to patient 11: **4** possibilities.
- Assign two of the three remaining N to patient 12: **3** possibilities; patient 13 receives the last N.

`4 × 15 × 4 × 3 = 720` valid N/C assignments consistent with the published table. Each assignment deterministically yields a donor-equal-weighted mean contrast `D=(D12+D13)/2`. Enumerating all assignments produces *identification bounds*:

| Cell population | Gene | Smallest possible donor-average (C − N) | Largest possible donor-average (C − N) | Direction invariant? |
|---|---|---:|---:|---|
| Basal keratinocytes | GATA6 | -0.05635 | -0.00785 | **Yes, lower in C** |
| Basal keratinocytes | FABP5 | +0.02139 | +0.21010 | **Yes, higher in C** |
| Basal keratinocytes | CRABP2 | -0.06220 | -0.01142 | **Yes, lower in C** |
| Transitional keratinocytes | GATA6 | -0.10775 | +0.01121 | **No** |
| Transitional keratinocytes | FABP5 | +0.39412 | +1.02888 | **Yes, higher in C** |
| Transitional keratinocytes | CRABP2 | -0.00940 | +0.39860 | **No** |
| HF keratinocytes | GATA6 | -0.36452 | +0.09878 | **No** |
| HF keratinocytes | FABP5 | -0.03959 | +0.89697 | **No** |

These are deterministic **extrema over assignments** using the deposited normalized/log expression values. They are NOT confidence intervals, p-values, fold changes, effect sizes on clinical acne, estimated probabilities of correctness, or proof of direction of causality. The 720 assignments need not be equally plausible biologically. Any verified partial sample-to-patient mapping would narrow these sets further. As this is only two paired donors for the C-N comparison and each donor's samples may come from different body sites, the inference is limited.

A particularly useful finding is that the **GATA6 decrease and FABP5 increase in basal keratinocytes** cannot be reversed merely by reassigning the *recent* N/C samples to patients in any manner consistent with Table S1, if the donor-averaging analysis is applied as defined. This is a robustness result, not independently validated mechanism.

### Exact executable analysis

Source: [round7_partial_donor_identification.py](round7_partial_donor_identification.py). Output: `gse301280_annotated_results/round7_donor_uncertainty_bounds.csv`. GitHub Actions additionally reruns it with internal assertions confirming 720 assignments and the three predefined invariant signs. Inputs are already-committed section-by-celltype summaries from the full h5ad.

This approach is an extension of partial-identification/bounds logic, **not a new general mathematical theorem**. It doesn't replace knowing the donor crosswalk, and conditions such as same donor/site matter for future causal interpretation.

## 3. An experimentally usable morphology model ALREADY exists

The apparently missing morphologic model is not wholly absent. **Guy, Green and Kealey (1996, J Invest Dermatol 106:176–182; doi:10.1111/1523-1747.ep12329907)** microdissected human sebaceous pilosebaceous infundibula and maintained them structurally for 7 days. They reported:

- controlled IL-1α stimulation caused **comedone-like hypercornification**;
- IL-1 receptor antagonist blocked the effect;
- EGF or TGF-α caused disorganization and apparent rupture of the isolated duct;
- retinoid exposure affected DNA synthesis/parakeratosis;
- about 20% of donors' isolated infundibula developed spontaneous hypercornification.

Their 1998 follow-up (doi:10.1159/000017862; PMID 9557221) also reports preserved viability and cytokine-associated hypercornification. A 1998 cytokine paper (PMID 9540984) distinguishes IL-1α-associated hypercornification from other cytokine-induced inflammatory surface markers.

**Strength:** This is a demonstrated human tissue model with a **structural** hypercornification endpoint, better for the first causal experiment than a sebocyte-only spheroid or fixed Xenium sections.

**Limitations:** It models histologic ductal changes in isolated human skin appendages over about **one week**. It does not fully reproduce human acne, variable skin microbiome, circulating endocrine signals, an intact lipid-providing sebaceous gland, years of comedogenesis or long-term off-treatment remission. Some donor tissues come from anatomy other than adolescent acne-prone face. Our lab has NOT recreated this model.

A 2026 primary-sebocyte sebonoid methods paper (doi:10.1016/j.jid.2026.04.007) is useful for **sebaceous target and lipid studies**, but not as a replacement for this infundibular plugging assay.

## 4. Precedent about remission that affects novelty

**Downie, Sanders and Kealey (2002, Br J Dermatol, doi:10.1046/j.1365-2133.2002.04946.x)** investigated human sebaceous gland responses to acne-relevant cytokines and described cellular changes that might help explain individual lesion remission. That is prior art for in-vitro modelling of remission and for cytokine-mediated plasticity.

Thus the idea 'acne lesions can remit via changes in sebocyte state' is not new. A novel, falsifiable FAR-2 contribution, if any, would be **specific, on-target reversibility of an acne-relevant *infundibular morphologic state* after verified withdrawal, distinguished from ordinary slow relaxation and reversible without injury**.

## 5. An incremental, decision-oriented experimental program

### Gate A — reproduce the classical model first (specialist academic lab)

Compare donor-matched isolated infundibula under vehicle, a validated comedogenic cytokine challenge and an established mechanistic blocking control, while monitoring blinded **infundibular hypercornification**, lumen plugging, viability and tissue integrity. Independently reproduce the historical morphology finding and measure inherent between-donor variability. A failed positive control means the model is not usable to evaluate FAR-2.

Avoid equating isolated ductal keratin accumulation with clinically meaningful, persistent human acne.

### Gate B — establish direct target necessity and rescue

Randomize viable donor-matched tissue units to target-specific GATA6 normalization versus appropriate non-targeting controls, and separately to FABP5 normalization versus controls. Include independent specificity-rescue approaches and a suitable standard acne-related retinoid comparator. Confirm protein-level perturbation and differentiated-cell function. **Pre-specify a primary structural occlusion endpoint**; analyze gene biomarkers secondarily.

Sufficient tissue quantity, delivery method, transduction artifacts, loss of cell viability and use of a multi-day model are critical feasibility barriers. Escalation of pleiotropic GATA6/TGFβ could disrupt barrier, immune regulation or tissue renewal. Direct FABP5 effects may vary with the lipid environment and retinoid receptor context; do not assume suppression always benefits the model.

**Decision:** Advance a target only if it changes morphology reproducibly in independent donors, the effect is target-specific and reverses with molecular rescue, and viability/barrier/inflammation safeguards pass. Otherwise do not claim a causal acne mechanism.

### Gate C — test interaction, timing and morphology

Once both single-target mechanisms can be manipulated cleanly, perform factorial experiments with consistent treatment timing relative to challenge. Compare pretreatment prevention, intervention after an established keratin plug, and matched washout. Track the relationship of IL-1α, TGFβ, lipid response and retinoid receptor activity with actual occlusion morphology.

A cytokine-driven positive control should **not** automatically be interpreted as faithfully recapitulating endocrine- or microbe-dependent acne initiation.

### Gate D — prove off-treatment durability only in a suitable longer-lived model

A 7-day organ culture might show *prevention or transient rescue* but cannot establish the kind of lasting follicular state reset FAR-2 envisions. For durability, require a viable, independent-donor human system that survives a sufficiently long drug-clearance/epithelial turnover/re-challenge period. Reversible perturbation, verified clearance, matched background inputs and slow-relaxation controls remain essential.

### Optional observational translation

In a separate ethically approved, non-interventional study, use serial LC-OCT or reflectance confocal microscopy to trace follicles from nonlesional to microcomedone to clinical acne. A 2019 longitudinal in-vivo study followed 10 people by RCM/D-OCT once weekly for six weeks (doi:10.1111/jdv.15641); this constitutes precedent, not novel first-in-human imaging. Such imaging could test prospective lesion-risk biomarkers, not validate untested compounds on people.

## 6. Decisions and remaining questions

- **Known source question solved partly:** We know the precise donor-level *group counts* from the actual supplement; we no longer need to invent a donor structure.
- **New computational advance:** Exact sign bounds for an equal-donor contrast under every category-compatible assignment show three consistent signals and several non-robust ones. This does not solve causal orientation.
- **Experimental model priority clarified:** Validate a human infundibular IL-1α hypercornification assay before elaborate sebonoid-only drug discovery or large virtual screens.
- **Real blockers:** an experimentally tractable on-target GATA6/FABP5 manipulation in viable isolated human ducts, sufficient matched tissue from independent donors, a reproducible quantitative occlusion endpoint, reliable safety measures and a longer-lived model for remission.
- **No human self-experimentation:** All biological interventions would be performed under qualified institutional research oversight. Our work here is reproducible computational evidence and a proposal, not a therapy.

## Primary source references

- Durgin et al. 2026 and Supplemental Table S1: https://insight.jci.org/articles/view/198021/sd/pdf/render/1
- Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280
- Guy, Green, Kealey 1996, human infundibular model: https://pubmed.ncbi.nlm.nih.gov/8592071/
- Guy and Kealey 1998 replication: https://pubmed.ncbi.nlm.nih.gov/9557221/
- Guy and Kealey 1998 cytokine specificity: https://pubmed.ncbi.nlm.nih.gov/9540984/
- Downie et al. 2002 human sebaceous remission model: https://pubmed.ncbi.nlm.nih.gov/12410695/
- Oulès et al. 2020 human follicle GATA6: https://www.nature.com/articles/s41467-020-18784-z
- Durgin et al. 2026 human AP-1/FABP5: https://insight.jci.org/articles/view/198021
- Sebonoids 2026 methods: https://pubmed.ncbi.nlm.nih.gov/42274454/
- Longitudinal RCM/OCT (2019): https://pubmed.ncbi.nlm.nih.gov/31025401/
