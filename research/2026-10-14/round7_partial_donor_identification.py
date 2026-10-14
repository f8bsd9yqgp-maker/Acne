"""Round 7: exact donor-mapping uncertainty bounds from original study's donor table.
PUBLIC cross-sectional research data; 720 FEASIBLE category assignments, NOT probabilities.

The dataset's *actual patient ID* is missing. We use the count constraints in the
published Supplement Table S1 (Durgin et al., JCI Insight 2026) WITHOUT INVENTING
a sample-to-patient crosswalk. The only patients with both recent nonlesional
and comedonal sections are patients 12 (2 N, 1 C) and 13 (1 N, 3 C).
Patients 10 (2 N, 1 P) and 11 (1 N, 1 P) supply the other 3 N.
Healthy donor sections and pustular sections are not used in this contrast.

For every assignment of N=6 and C=4 relevant sections compatible with those
counts, we calculate, for each gene and cell type, the EQUAL-DONOR average of
the within-patient mean C-N contrast across donors 12 and 13.

An invariant sign across 720 arbitrary category-compatible assignments DOES
NOT establish statistical significance, longitudinal order, causal direction,
sampling-site comparability, or that each possible assignment is equally likely.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
BASE=ROOT/"gse301280_annotated_results"
DATA=BASE/"annotated_sample_celltype_gene_expression.csv"
OUT=BASE/"round7_donor_uncertainty_bounds.csv"
CATS={"non":"Acne Non-lesional","com":"Acne Comedo"}
GENES=("GATA6","FABP5","CRABP2","KRT79","KRT10","RARA","RXRA")
TYPES=("Basal KC","Transitional KC","HF KC","Spinous KC","Sebocyte")
def enumerate_compatible(non,com):
    """p10:2N, p11:1N, p12:2N+1C, p13:1N+3C."""
    non=tuple(sorted(non)); com=tuple(sorted(com))
    assert len(non)==6 and len(com)==4,(non,com)
    for c12 in com:
        c13=tuple(s for s in com if s != c12)
        for n10 in combinations(non,2):
            rest1=tuple(s for s in non if s not in n10)
            for n11 in rest1:
                rest2=tuple(s for s in rest1 if s != n11)
                for n12 in combinations(rest2,2):
                    n13=tuple(s for s in rest2 if s not in n12)
                    assert len(n13)==1
                    yield {"12_N":n12, "12_C":(c12,), "13_N":n13,"13_C":c13}
def main():
    df=pd.read_csv(DATA)
    req={"sample","disease","is_recent_sample","cell_type","gene","n_cells","mean_anndata_X"}
    assert req.issubset(df.columns)
    recent=df[df.is_recent_sample.astype(str).str.lower()=="true"]
    unique=recent[["sample","disease"]].drop_duplicates()
    non=tuple(unique.loc[unique.disease==CATS["non"],"sample"])
    com=tuple(unique.loc[unique.disease==CATS["com"],"sample"])
    assert len(non)==6 and len(com)==4
    assignments=list(enumerate_compatible(non,com))
    assert len(assignments)==720
    result=[]
    for typ in TYPES:
        for gene in GENES:
            q=recent[(recent.cell_type==typ)&(recent.gene==gene)]
            assert not q.sample.duplicated().any(),(typ,gene)
            sig=dict(zip(q["sample"],q["mean_anndata_X"]))
            cov=dict(zip(q["sample"],q["n_cells"]))
            if any(v not in sig for v in non+com):continue
            vals12=[];vals13=[];av=[]
            for a in assignments:
                mu=lambda keys:sum(sig[key] for key in keys)/len(keys)
                c12=mu(a["12_C"])-mu(a["12_N"])
                c13=mu(a["13_C"])-mu(a["13_N"])
                vals12.append(c12);vals13.append(c13);av.append((c12+c13)/2)
            av=np.array(av);vals12=np.array(vals12);vals13=np.array(vals13)
            if (av>1e-12).all(): verdict="POSITIVE_FOR_ALL_COMPATIBLE_ASSIGNMENTS"
            elif (av< -1e-12).all():verdict="NEGATIVE_FOR_ALL_COMPATIBLE_ASSIGNMENTS"
            else:verdict="AMBIGUOUS"
            result.append(dict(
                cell_type=typ,gene=gene,eligible_nonlesional_sections=6,
                eligible_comedonal_sections=4,compatible_assignments=len(assignments),
                donor12_delta_min=float(vals12.min()),donor12_delta_max=float(vals12.max()),
                donor13_delta_min=float(vals13.min()),donor13_delta_max=float(vals13.max()),
                equal_donor_average_min=float(av.min()),
                equal_donor_average_max=float(av.max()),
                positive_assignments=int(np.count_nonzero(av>0)),
                negative_assignments=int(np.count_nonzero(av<0)),
                verdict=verdict,
                minimum_cells_any_section=min(cov[v] for v in non+com),
                limitations="Bounds over category-compatible mappings; not probability/significance/causality."
            ))
    o=pd.DataFrame(result)
    o.to_csv(OUT,index=False)
    selected=o[o.gene.isin(["GATA6","FABP5","CRABP2"]) &
        o.cell_type.isin(["Basal KC","Transitional KC","HF KC"])]
    print("Exhaustive category-compatible assignment count: 720")
    print(selected[["cell_type","gene","equal_donor_average_min",
      "equal_donor_average_max","verdict","minimum_cells_any_section"]].to_string(index=False))
    for typ,g,sgn in [("Basal KC","GATA6","NEGATIVE_FOR_ALL_COMPATIBLE_ASSIGNMENTS"),
                      ("Basal KC","FABP5","POSITIVE_FOR_ALL_COMPATIBLE_ASSIGNMENTS"),
                      ("Transitional KC","FABP5","POSITIVE_FOR_ALL_COMPATIBLE_ASSIGNMENTS")]:
        row=o[(o.cell_type==typ)&(o.gene==g)]
        assert len(row)==1 and row.iloc[0].verdict==sgn,(typ,g,row)
    print("IMPORTANT: Actual patient mapping remains unresolved. Sign does not equal causation.")
if __name__=="__main__":main()
