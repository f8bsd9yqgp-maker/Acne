"""Exact partial identification across ALL THREE donor-matched N/C patients.

Source: Durgin et al 2026 Supplement Table S1.
- Patient 7 archival: 1 N, 3 C.
- Patient 12 recent: 2 N, 1 C.
- Patient 13 recent: 1 N, 3 C.
- Patient 8 + patient 9 archival: each 1 N, no C.
- Patients 10/11 recent: 2 N + 1 N respectively, no C.

Actual sample↔patient matching remains unknown. We enumerate 3 options for
patient 7's 1 N section multiplied by 720 options for 12/13 and the two
unpaired recent donors = 2,160 assignments. Contrast is equal-weighted
average of the 3 paired donor-level C-minus-N mean differences.

The bounds are exact for THIS descriptive statistic under category constraints,
not a confidence interval, posterior probability, significance test, or
evidence about causation or lesion persistence. Biopsy sites differ.
"""
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from round7_partial_donor_identification import enumerate_compatible

ROOT=Path(__file__).resolve().parent
BASE=ROOT/"gse301280_annotated_results"
GENES=("GATA6","FABP5","CRABP2","KRT5","KRT79","KRT10")
CELLTYPES=("Basal KC","Transitional KC","HF KC","Spinous KC","Sebocyte")
def main():
    dat=pd.read_csv(BASE/"annotated_sample_celltype_gene_expression.csv")
    samplemeta=(dat.loc[dat.gene=="GATA6",["sample","disease","is_recent_sample"]]
                .drop_duplicates())
    get=lambda recent,cat:tuple(samplemeta.loc[
        (samplemeta.is_recent_sample==recent)&(samplemeta.disease==cat),"sample"])
    arch_N=get(False,"Acne Non-lesional")
    arch_C=get(False,"Acne Comedo")
    rec_N=get(True,"Acne Non-lesional")
    rec_C=get(True,"Acne Comedo")
    assert tuple(map(len,(arch_N,arch_C,rec_N,rec_C)))==(3,3,6,4)
    rec_assign=list(enumerate_compatible(rec_N,rec_C))
    assert len(rec_assign)==720
    totals=[]
    for cell_type in CELLTYPES:
        for gene in GENES:
            x=dat[(dat.cell_type==cell_type)&(dat.gene==gene)]
            vals=dict(zip(x["sample"],x["mean_anndata_X"]))
            count=dict(zip(x["sample"],x["n_cells"]))
            needed=arch_N+arch_C+rec_N+rec_C
            assert all(k in vals for k in needed)
            mu=lambda ids: sum(vals[k] for k in ids)/len(ids)
            d7=[mu(arch_C)-mu([a]) for a in arch_N]
            d12=[];d13=[]
            for a in rec_assign:
                d12.append(mu(a["12_C"])-mu(a["12_N"]))
                d13.append(mu(a["13_C"])-mu(a["13_N"]))
            # Exhaustive 2160 three-donor mean contrasts.
            full=np.array([(x+y+z)/3 for x in d7 for y,z in zip(d12,d13)])
            assert full.shape==(2160,)
            if (full>0).all():sign="ALWAYS_POSITIVE"
            elif (full<0).all():sign="ALWAYS_NEGATIVE"
            else:sign="AMBIGUOUS"
            totals.append(dict(cell_type=cell_type,gene=gene,
               archived_donor7_delta_min=min(d7),archived_donor7_delta_max=max(d7),
               donor_equal_C_minus_N_lower_bound=float(full.min()),
               donor_equal_C_minus_N_upper_bound=float(full.max()),
               possible_mappings=2160,sign_invariant=sign,
               n_archival_sections=6,n_recent_sections=10,
               min_annotated_cells_any_section=min(count[k] for k in needed)))
    result=pd.DataFrame(totals)
    result.to_csv(BASE/"round7_three_donor_partial_identification.csv",index=False)
    show=result[result.gene.isin(["GATA6","FABP5","CRABP2"])&
         result.cell_type.isin(["Basal KC","Transitional KC","HF KC"])]
    print(show[["cell_type","gene","donor_equal_C_minus_N_lower_bound",
           "donor_equal_C_minus_N_upper_bound","sign_invariant",
           "min_annotated_cells_any_section"]].to_string(index=False))
    def observed(ct,g):return result[(result.cell_type==ct)&(result.gene==g)].iloc[0]
    for ct in ["Basal KC","Transitional KC"]:
        assert observed(ct,"FABP5").sign_invariant=="ALWAYS_POSITIVE"
    assert observed("Basal KC","GATA6").sign_invariant=="AMBIGUOUS"
    print("VERIFIED all 2160 category-compatible mappings. No actual donor crosswalk inferred.")
if __name__=="__main__":main()
