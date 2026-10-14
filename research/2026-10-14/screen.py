"""Hypothesis-generator for acne intervention targets, NOT molecule screening or a validated human model.
Reproducible deterministic grid screen. Standard Python, NumPy, Pandas, Matplotlib.
"""
from itertools import combinations, product
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
TARGETS = ["lipid_quality", "ductal_identity", "GehA_neutralization", "strain_ecology", "inflammation_resolution", "barrier_support", "sebum_normalization", "retinoid_routing"]
PHENOTYPES = {
    "comedonal":   dict(seb=.74, duct=.85, micro=.40, inflam=.30, lipid=.68, barrier=.65),
    "inflammatory":dict(seb=.67, duct=.60, micro=.82, inflam=.85, lipid=.74, barrier=.54),
    "oily_sensitive":dict(seb=.92, duct=.62, micro=.64, inflam=.65, lipid=.88, barrier=.26),
    "fragile_barrier":dict(seb=.57, duct=.53, micro=.58, inflam=.72, lipid=.72, barrier=.13),
}
# Explicit theoretical uncertainty around mechanistic interaction weights.
WORLDS = [
    dict(lipid=1, duct=1, geh=1, immune=1),
    dict(lipid=.7, duct=1.3, geh=.6, immune=1.3),
    dict(lipid=1.3, duct=.7, geh=1.5, immune=.7),
    dict(lipid=.85, duct=1.15, geh=1.25, immune=1.15),
]

def candidates():
    # One-, two-, and three-target combinations at 10 nonzero control levels;
    # four-target combinations at 5 nonzero control levels.
    batches = [np.zeros((1, len(TARGETS)), dtype=np.float32)]
    for k in (1, 2, 3, 4):
        levels = np.linspace(.1, 1.0, 10, dtype=np.float32) if k < 4 else np.linspace(.2, 1.0, 5, dtype=np.float32)
        levels_grid = np.array(list(product(levels, repeat=k)), dtype=np.float32)
        for indices in combinations(range(len(TARGETS)), k):
            a = np.zeros((len(levels_grid), len(TARGETS)), dtype=np.float32)
            a[:, indices] = levels_grid
            batches.append(a)
    return np.concatenate(batches)

def evaluate(controls, phenotype, world):
    """Phenomenological fixed-point risk proxy, no kinetic rates or patient calibration."""
    L, D, G, M, I, B, S, R = [controls[:, j] for j in range(8)]
    p = PHENOTYPES[phenotype]
    w = world
    # Basal physiological sebum always retained; too much suppression harms barrier.
    seb = (.20 + .80*p['seb'])*(1-.46*S)
    duct = np.full_like(S, p['duct'])
    infl = np.full_like(S, p['inflam'])
    lipid = np.full_like(S, p['lipid'])
    micro = np.full_like(S, p['micro'])
    for _ in range(16):
        new_lipid = np.clip((.55*p['lipid'] + w['lipid']*(.22*seb + .13*infl))*(1-.60*L),0,2)
        new_micro = np.clip((.68*p['micro'] + .20*new_lipid + .12*seb)*(1-.64*M),0,2)
        geh = np.clip(w['geh']*.76*new_micro*(.40 + .60*new_lipid)*(1-.80*G), 0, 2)
        new_duct = np.clip((.55*p['duct'] + w['duct']*(.22*new_lipid+.20*infl+.11*seb))*(1-.53*D)*(1-.20*R),0,2)
        new_infl = np.clip((.35*p['inflam']+w['immune']*(.23*geh+.22*new_lipid+.22*new_duct+.11*new_micro))*(1-.50*I)*(1-.11*B),0,2)
        duct = .40*duct+.60*new_duct
        infl = .40*infl+.60*new_infl
        lipid = .40*lipid+.60*new_lipid
        micro = .40*micro+.60*new_micro
    # Composite is a *research heuristic*, not a calibrated acne or clinical efficacy score.
    risk = .40*duct + .36*infl + .14*lipid + .10*geh
    barrier_integrity = np.clip(p['barrier'] + .34*B - .14*infl - .13*lipid -.17*np.maximum(S-.55,0) - .11*np.maximum(R-.55,0),0,1)
    # Safety uncertainty proxy; values and weights entirely assumed.
    toxicity = (.090*S**2 + .060*I**2 + .050*R**2 + .033*D**2 + .018*M**2 + .010*G**2 + .009*L**2 + .004*B**2)
    return risk.astype(np.float32), barrier_integrity.astype(np.float32), toxicity.astype(np.float32)

def main():
    c = candidates()
    n = c.shape[0]
    print('candidate_count',n)
    assert n == 102631, (n, 'Unexpected grid count')
    # Compare 16 phenotype × uncertainty scenarios against untreated reference.
    reductions, barriers, toxs = [], [], []
    control = np.zeros((1, len(TARGETS)), dtype=np.float32)
    for phenotype in PHENOTYPES:
        for world in WORLDS:
            risk, barrier, tox = evaluate(c, phenotype, world)
            b0 = evaluate(control, phenotype, world)[0][0]
            reductions.append((b0-risk)/max(b0, 1e-6))
            barriers.append(barrier)
            toxs.append(tox)
    reductions = np.array(reductions)
    barriers = np.array(barriers)
    toxs = np.array(toxs)
    mean_improvement = reductions.mean(axis=0)
    worst_improvement = reductions.min(axis=0)
    # Barrier health must not deteriorate more than 0.03 from a same-phenotype same-world untreated reference.
    barrier_changes = []
    for phenotype in PHENOTYPES:
        for world in WORLDS:
            base_barrier = evaluate(control, phenotype, world)[1][0]
            change = barriers[len(barrier_changes)] - base_barrier
            barrier_changes.append(change)
    worst_barrier_change = np.array(barrier_changes).min(axis=0)
    mean_toxicity=toxs.mean(axis=0)
    # Relative heuristic: rewarded robust improvement, penalized toxicity, plus simplicity.
    num_targets = (c > 1e-6).sum(axis=1)
    score = .55*worst_improvement+.45*mean_improvement -1.15*mean_toxicity - .009*np.maximum(num_targets-1,0)
    safe = (worst_barrier_change>=-.03)&(mean_toxicity<=.055)
    score[~safe] = -999
    df = pd.DataFrame(c,columns=TARGETS)
    df['number_of_targets']=num_targets
    df['mean_proxy_reduction']=mean_improvement
    df['worst_proxy_reduction']=worst_improvement
    df['min_barrier_change']=worst_barrier_change
    df['toxicity_proxy']=mean_toxicity
    df['screen_score']=score
    df['passes_heuristic_guardrails']=safe
    df.to_csv(OUT/'all_virtual_interventions.csv',index=False)
    df.sort_values('screen_score',ascending=False).head(40).to_csv(OUT/'top_40_virtual_interventions.csv',index=False)
    # Package best candidate family per distinct set of active targets, not duplicates with slightly different strengths.
    passed = df.loc[safe].copy()
    passed['family']=[' + '.join([TARGETS[j] for j in np.flatnonzero(row)]) or 'untreated' for row in c[safe]]
    family = passed.sort_values('screen_score',ascending=False).drop_duplicates('family')
    family.head(25).to_csv(OUT/'top_25_distinct_target_combinations.csv',index=False)
    # Single target high-level signal and flags.
    single = df.loc[df['number_of_targets']==1]
    single_best = single.sort_values('screen_score',ascending=False).drop_duplicates(TARGETS).head(20)
    single_best.to_csv(OUT/'single_target_comparators.csv',index=False)
    summary={
        'screened_interventions':n, 'safe_interventions':int(safe.sum()),
        'safe_share':float(safe.mean()), 'best_score':float(family.iloc[0]['screen_score']),
        'best_mean_reduction':float(family.iloc[0]['mean_proxy_reduction']),
        'best_worst_reduction':float(family.iloc[0]['worst_proxy_reduction']),
        'best_family':family.iloc[0]['family'],
    }
    (OUT/'screen_summary.txt').write_text('\n'.join(f'{k}: {v}' for k,v in summary.items())+'\n')
    print(summary)
    print('\nLeading distinct strategies:')
    print(family[['family','mean_proxy_reduction','worst_proxy_reduction','min_barrier_change','toxicity_proxy','screen_score']].head(12).to_string(index=False))
    print('\nBest individual-target examples:')
    print(single.sort_values('screen_score',ascending=False).assign(family=lambda d:[TARGETS[np.flatnonzero(c[i])[0]] for i in d.index]).drop_duplicates('family')[['family','mean_proxy_reduction','screen_score']].to_string(index=False))
    # Charts without fixed coloring
    chart=family.head(10).iloc[::-1]
    plt.figure(figsize=(10,5.7))
    plt.barh(np.arange(len(chart)),chart['worst_proxy_reduction']*100,alpha=.9,label='Worst of 16 scenarios')
    plt.yticks(np.arange(len(chart)),chart['family'],fontsize=8)
    plt.xlabel('Simulated acne-risk proxy reduction (%) — NOT clinical efficacy')
    plt.title('Top distinct intervention combinations: assumed model')
    plt.tight_layout()
    plt.savefig(OUT/'intervention_ranking.png',dpi=165)
    plt.close()
    plt.figure(figsize=(7.5,5.5))
    sub=df.loc[safe].sample(n=min(5000,int(safe.sum())),random_state=2026)
    plt.scatter(sub['toxicity_proxy'],sub['mean_proxy_reduction']*100,alpha=.17,s=11)
    plt.scatter(family.iloc[0]['toxicity_proxy'],family.iloc[0]['mean_proxy_reduction']*100,marker='x',s=100,linewidths=2)
    plt.xlabel('Assumed safety cost (arbitrary units; lower better)')
    plt.ylabel('Simulated mean risk-proxy reduction (%) — NOT measured')
    plt.title('Illustrative efficacy–safety trade-off')
    plt.tight_layout()
    plt.savefig(OUT/'safety_tradeoff.png',dpi=165)
    plt.close()
    return summary

if __name__=='__main__':main()
