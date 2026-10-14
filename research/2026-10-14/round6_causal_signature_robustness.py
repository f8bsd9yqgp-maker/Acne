"""Round 6: sample-level, descriptive robustness of target-gene patterns.

Uses already derived GSE301280 annotated per-sample per-celltype measurements.
This cannot identify causal direction or infer donor-independent significance.
No p-values: sections from the same person may be correlated.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
SRC=ROOT/"gse301280_annotated_results"/"annotated_sample_celltype_gene_expression.csv"
OUT=ROOT/"gse301280_annotated_results"/"round6_section_level_robustness.csv"
CATEGORIES=["Healthy","Acne Non-lesional","Acne Comedo","Acne Pustule"]
CONTRASTS=[("Healthy","Acne Non-lesional","acne_prone_nonlesional_vs_healthy"),
 ("Acne Non-lesional","Acne Comedo","comedone_vs_acne_prone_nonlesional"),
 ("Acne Non-lesional","Acne Pustule","pustule_vs_acne_prone_nonlesional")]
TARGETS=["GATA6","FABP5","CRABP2","KRT79","KRT5","KRT10","RARA","RXRA"]
CELLTYPES=["Basal KC","Transitional KC","HF KC","Sebocyte","Spinous KC"]

def median(x):
    return float(np.median(x)) if len(x) else np.nan

def cliffs_delta(x,y):
    if not len(x) or not len(y):return np.nan
    x=np.asarray(x);y=np.asarray(y)
    return float((np.sum(x[:,None]>y[None,:])-np.sum(x[:,None]<y[None,:]))/(len(x)*len(y)))

def loo(x,y):
    if min(len(x),len(y))<2:return (np.nan,np.nan,False)
    differences=[median(np.delete(x,i))-median(y) for i in range(len(x))]
    differences +=[median(x)-median(np.delete(y,j)) for j in range(len(y))]
    observed=median(x)-median(y)
    # Median contrast sign robustness is a technical sensitivity, not significance.
    stable=(observed!=0) and all(np.sign(d)==np.sign(observed) for d in differences)
    return (float(min(differences)),float(max(differences)),bool(stable))

def main():
    df=pd.read_csv(SRC)
    assert {"sample","disease","is_recent_sample","cell_type","gene","mean_anndata_X","n_cells"}.issubset(df.columns)
    meta=df.loc[df.gene=="GATA6",["sample","disease","is_recent_sample"]].drop_duplicates()
    assert meta.groupby("disease")["sample"].nunique().reindex(CATEGORIES).tolist()==[5,9,7,3]
    out=[]
    for cell in CELLTYPES:
      for gene in TARGETS:
        sub=df[(df.cell_type==cell)&(df.gene==gene)]
        for filter_name in ("all", "recent", "recent_min100"):
          chosen=sub.copy()
          if "recent" in filter_name:
             chosen=chosen[chosen.is_recent_sample.astype(str).str.lower()=="true"]
          if "min100" in filter_name:
             chosen=chosen[chosen.n_cells>=100]
          for base,case,label in CONTRASTS:
            y=chosen.loc[chosen.disease==base,"mean_anndata_X"].to_numpy(dtype=float)
            x=chosen.loc[chosen.disease==case,"mean_anndata_X"].to_numpy(dtype=float)
            if len(x)==0 or len(y)==0:continue
            lo,hi,stable=loo(x,y)
            out.append({"cell_type":cell,"gene":gene,"contrast":label,
                        "subset":filter_name,"n_case_sections":len(x),"n_reference_sections":len(y),
                        "case_median":median(x),"reference_median":median(y),
                        "delta_median_X":median(x)-median(y),"cliffs_delta":cliffs_delta(x,y),
                        "loo_min_delta":lo,"loo_max_delta":hi,
                        "loo_sign_stable":stable,
                        "interpretation":"descriptive sections; not independent donors, not causation"})
    result=pd.DataFrame(out)
    result.to_csv(OUT,index=False)
    show=result[(result.cell_type.isin(["Basal KC","Transitional KC","HF KC"])) &
        (result.gene.isin(["GATA6","FABP5"])) &
        (result.contrast=="comedone_vs_acne_prone_nonlesional") &
        (result.subset=="recent_min100")]
    print(show[["cell_type","gene","n_case_sections","n_reference_sections",
       "case_median","reference_median","delta_median_X","cliffs_delta",
       "loo_min_delta","loo_max_delta","loo_sign_stable"]].to_string(index=False))
    print("\nCAUTION: Same-patient samples and batch effects invalidate section-independent statistical inference.")
if __name__=="__main__":main()
