"""Alternative hypothetical objective functions. These are assumptions, NOT validations."""
import pandas as pd
from pathlib import Path
P=Path(__file__).parent
D=pd.read_csv(P/'all_virtual_interventions.csv')
ctrl=['lipid_quality','ductal_identity','GehA_neutralization','strain_ecology','inflammation_resolution','barrier_support','sebum_normalization','retinoid_routing']
safe=D['passes_heuristic_guardrails'].astype(bool)
D=D[safe].copy()
D['family']=D[ctrl].gt(0).apply(lambda r:' + '.join(r.index[r]) if r.any() else 'untreated',axis=1)
scenario_scores={
    'baseline':D['screen_score'],
    'safety_first':.50*D['mean_proxy_reduction']+.50*D['worst_proxy_reduction']-2.0*D['toxicity_proxy']+.10*D['min_barrier_change']-.015*(D['number_of_targets']-1).clip(lower=0),
    'simplicity_first':.55*D['worst_proxy_reduction']+.45*D['mean_proxy_reduction']-1.15*D['toxicity_proxy']-.06*(D['number_of_targets']-1).clip(lower=0),
    'immune_caution':.55*D['worst_proxy_reduction']+.45*D['mean_proxy_reduction']-1.15*D['toxicity_proxy']-.009*(D['number_of_targets']-1).clip(lower=0)-.30*D['inflammation_resolution']**2,
}
res=[]
for label,score in scenario_scores.items():
    sub=D.copy();sub['alt_score']=score
    first=sub.sort_values('alt_score',ascending=False).drop_duplicates('family').iloc[0]
    res.append({'model_objective':label,'top_family':first['family'],'top_score_arbitrary':first['alt_score'],'mean_proxy_reduction':first['mean_proxy_reduction'],'worst_proxy_reduction':first['worst_proxy_reduction']})
answer=pd.DataFrame(res)
answer.to_csv(P/'objective_sensitivity.csv',index=False)
print(answer.to_string(index=False))
