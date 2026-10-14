# Round 3 — Follicular state switching, identifiability and model feasibility
**Research date:** 2026-10-14  
**Status:** New theoretical analysis + reproducible toy-model calculations, external literature check and research design. No new human data, lab experiment, validated drug or cure.

## A. Two obstacles that change our direction

**Obstacle 1 — phenotype validity:** We cannot conclude that we prevented comedones using only a sebaceous-gland organoid, an epidermal cell culture, or an isolated bacterial-growth assay. Each captures a different compartment. The desired *primary* research phenotype must be validated upper-infundibular keratinization/occlusion/microcomedone formation plus preservation of tissue function.

- Oulès et al. (2020) specifically performed GATA6 knockdown in **human hair follicle explants**. The authors cautioned that differentiation changes were modest because the knockdown kinetics and culture duration were too short for hyperkeratinization to develop: https://www.nature.com/articles/s41467-020-18784-z
- A June 2026 sebaceous organoid methods report describes primary-human-derived **sebonoids** useful for gland differentiation, but it does **not** establish the upper-infundi­bular comedogenesis endpoint: https://www.sciencedirect.com/science/article/pii/S0022202X26010511
- A 2025 **SCmic** microecology system can evaluate multi-strain bacterial communities and lipids but deliberately uses **nonviable** epidermal/sebaceous cells. It therefore cannot demonstrate living keratinocyte switching or regenerative durability: https://www.nature.com/articles/s41467-025-58377-2
- Durgin et al. (2026) observed **different** disease-dependent KRT10 responses: elevated in comedones, often lower in pustules. A single scalar 'normalization' metric that treats all lesions alike may be biologically misleading: https://insight.jci.org/articles/view/198021

**Obstacle 2 — dynamic identifiability:** An apparently persistent treatment effect can be (i) residual pharmacology, (ii) slow unique-equilibrium recovery, (iii) cell population selection, (iv) long-lived niche signaling or (v) true multistability. Our old 102,631-profile fixed-point screen has **no time, differentiation/turnover history or drug clearance dynamics**; it cannot distinguish any of these. Its rank scores must not be used to claim a reset.

## B. Exact mathematical result in a minimal illustrative model

Consider the normalized canonical one-dimensional cubic:
`dx/dt = x - x^3 + h(t)`.

**This is a textbook bistability illustration, not a derived equation for real human follicles.**

For constant h, fixed points satisfy `h=x^3-x`, and stability requires `1-3x^2<0`. There are two stable branches when:
`|h|<h_c=2/(3 sqrt(3)) ≈ 0.3849001794597505`.

If input h is swept quasistatically across both fold bifurcations, the limiting hysteresis loop area is **exactly 3/2** in arbitrary x–h units. To see this, let `F(x)=(3/4)x^4-(1/2)x^2`, so `∫x dh = F(x)` along equilibrium branches. Each branch integral between fold limits has magnitude 3/4; their difference is 3/2.

Compare with a single-attractor lag model:
`dx/dt=(h(t)-x)/tau`, where tau>0. At a fixed input h=0, differences in initial states shrink exactly as `Δx(t)=Δx(0)exp(-t/tau)`. However, if tau is very large relative to the experiment duration, finite-duration traces imitate persistent state separation. For any finite observation T and fractional tolerance epsilon in (0,1), choosing `tau> -T/ln(1-epsilon)` makes the contraction fraction **smaller** than epsilon. Thus no finite waiting-time test can logically rule out *all* slow monostable alternatives without independent bounds.

The discriminating research strategy is **time-scale triangulation**, not a magic single endpoint:
- Verify active target engagement and the drug are absent.
- Repeatedly remeasure morphology and function under matched inputs.
- Compare hysteresis loops as external challenges ramp **increasingly slowly** while maintaining comparable tissue viability.
- Separately vary pulse duration and amplitude, track cell-identity shifts and cell-population turnover, and measure whether apparent states continue to diverge from background adaptation.
- Check two-way switching and recovery, not only 'after treatment looks better'.
- Even persistent loops *favor* rather than logically prove bistability because high-dimensional slow adaptation, altered cell composition or irreversibility can also create history dependence.

## C. Computational check (executed, deterministic)

Source: `state_switch_diagnostic.py`. Numerical ODE solver uses scipy `solve_ivp` with tight tolerances and max step 1. Synthetic input ramp goes -0.6 → +0.6 → -0.6, with tau=15 for single-attractor model. Duration is **per leg** and dimensionless.

| Ramp duration | Single-attractor loop area | Bistable cubic loop area |
|---|---:|---:|
| 20 | 0.527772 | 2.190261 |
| 80 | 0.389099 | 1.903746 |
| 320 | 0.125506 | 1.661771 |
| 1280 | 0.033154 | 1.563970 |

The cubic model approaches the derived asymptote of 1.5; the monostable loop contracts toward 0 as the ramp slows. These are **arbitrary simulation units**, not biomarkers, lesion counts, treatment response, human data, or an estimated experiment duration. Neither model was fitted or validated against acne observations.

`state_switch_diagnostic.py` also writes `washout_separation.csv` comparing theoretical state separation with tau=100 versus 1000. At time 100, the monostable separation is 0.73576 with tau=100, but 1.80967 with tau=1000, compared with a bistable value of 2. A single finite assay cannot identify an attractor.

## D. Concrete experimental roadmap

**Gate 0 — model eligibility:** Establish if a patient-derived follicular explant or reconstructed complete pilosebaceous system can reproducibly undergo **new upper-follicle plugging** while remaining viable long enough for a verified drug washout and delayed re-challenge. A short explant cannot answer an infinite-duration claim. A sebonoid or microbiome microecology model may supply *secondary* mechanistic assays, not replace primary morphology.

**Gate 1 — target and phenotype causality:** Compare independent orthogonal GATA6-related and FABP5/CRABP2-related perturbations with specificity rescue and reference retinoid in matched donor follicles. Normalized target gene levels are not sufficient without a functional structural endpoint; screen injury and barrier changes.

**Gate 2 — treatment-history challenge:** Randomize follicles within donor to preparation histories, confirm clearance, preserve identical endocrine/microbial background and record morphology at multiple durations. Separate comedonal and pustular phenotypes and the nonlesional high-risk population.

**Gate 3 — test for multistability:** At the same external conditions, test whether two paths produce reproducible stable endpoints, switching in both directions, pulse threshold behavior and rate-stable loop areas. Assess population selection via progenitor-marker and lineage-compatible assays. Endpoint is **model-specific history dependence**, not human remission.

**Gate 4 — real-world relevance:** Only after replicated causal effects and safety, develop testable delivery technology and contact independent investigators about human clinical relevance.

**Stop conditions:** if no suitable comedogenesis model exists; if effect only changes biomarkers; if pharmacology persists; if viability/barrier decreases; if the 'two states' disappear with adequate equilibration or cell-turnover controls. Do not call any of these a cure.

## E. Public human dataset identified but not reanalyzed
Durgin et al. deposited **GSE301280** at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280. It includes a ~758 MB combined `h5ad` file, with **24 tissue sections from 11 donors** and a targeted **100-gene panel**. The paper's original computational pipeline is https://github.com/joe-durgin/TranscriptFocusedSegmentation (reported commit `35b6c02`). This is a **cross-sectional** dataset; it cannot directly test future lesion incidence or response after drug withdrawal. A donor-aware, cell-compartment-aware replication of GATA6/KRT79, FABP5/CRABP2 and RAR targets is worthwhile. Technical cells from the same donor cannot substitute for donor replication. The dataset was *located* but its original h5ad file has **not** been downloaded or reprocessed here.

### Data analysis safeguards
- Check whether donor labels and sample age/archival status can be recovered before any donor-level claims.
- Stratify comedonal and pustular lesions instead of pooling.
- Use donor-level aggregation or defensible clustered/hierarchical inference, with paired sections analyzed as such.
- Ensure comparisons of FABP5/CRABP2 do not produce artifacts when CRABP2 is close to zero; compare expression by compartment and report detection frequencies.
- Confirm protein/pathway function independently of targeted 100-gene mRNA panels.

## F. Reproducibility bug identified and fixed
The `single_target_comparators.csv` export in the original archive incorrectly called `drop_duplicates(TARGETS)`, keeping the top 20 *doses* rather than the best result from each of 8 mechanisms; the old output only included **3 target types** (8 lipid quality, 7 ductal identity and 5 inflammation resolution rows). The GitHub `screen.py` script now assigns each active single-axis profile a `target` and deduplicates by target, asserting exactly eight comparator rows. The original ZIP remains unchanged as a historical artifact.

## Interpretation
**New original work in this update:** explicit identifiability comparison, reproducible numerical ramp calculations, exact canonical-model hysteresis-area derivation, model-eligibility evidence synthesis and a fixed export bug. **Not original biological proof:** generic bistability equations are established mathematical concepts. **Not obtained:** actual human molecular perturbation data, new compounds, a working full acne follicle model or an acne cure.
