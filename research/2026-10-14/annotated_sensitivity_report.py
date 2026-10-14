"""Postprocess precomputed public GSE301280 annotated per-section outputs.
No donor independence, causal interpretation, or human clinical inference.
Does NOT require redownloading the 758MB h5ad.
Run: python annotated_sensitivity_report.py
"""
from pathlib import Path
import pandas as pd

BASE=Path(__file__).resolve().parent/"gse301280_annotated_results"
GENES=["GATA6","FABP5","CRABP2"]
TYPES=["HF KC","Basal KC","Transitional KC","Sebocyte","Spinous KC"]
CATEGORIES=["Healthy","Acne Non-lesional","Acne Comedo","Acne Pustule"]
def main():
    df=pd.read_csv(BASE/"annotated_sample_celltype_gene_expression.csv")
    required={"sample","disease","cell_type","gene","n_cells","mean_anndata_X","is_recent_sample"}
    assert required.issubset(df.columns),required-set(df.columns)
    assert df["sample"].nunique()==24
    assert len(df[df.gene=="GATA6"].groupby(["sample","cell_type"]))==len(df[df.gene=="FABP5"].groupby(["sample","cell_type"]))
    assert df.groupby("sample")["disease"].nunique().max()==1
    metadata=df[["sample","disease","is_recent_sample"]].drop_duplicates()
    assert (metadata.groupby("disease")["sample"].nunique()
            .reindex(CATEGORIES).tolist()==[5,9,7,3])
    rows=[]
    for cell in TYPES:
        for gene in GENES:
            for disease in CATEGORIES:
                selected=df[(df.cell_type==cell)&(df.gene==gene)&(df.disease==disease)]
                for subset in ("all","recent","min100","recent_min100"):
                    sub=selected.copy()
                    if "recent" in subset:
                        sub=sub[sub.is_recent_sample.astype(str).str.lower()=="true"]
                    if "min100" in subset:sub=sub[sub.n_cells>=100]
                    rows.append({"cell_type":cell,"gene":gene,"disease":disease,
                        "filter":subset,"n_sections":len(sub),
                        "median_section_mean_expression":(sub.mean_anndata_X.median() if len(sub) else float("nan")),
                        "median_section_detection":(sub.fraction_nonzero.median() if len(sub) else float("nan")),
                        "min_cells_per_section":(int(sub.n_cells.min()) if len(sub) else None),
                        "median_cells_per_section":(float(sub.n_cells.median()) if len(sub) else None),
                        "max_cells_per_section":(int(sub.n_cells.max()) if len(sub) else None)})
    result=pd.DataFrame(rows)
    result.to_csv(BASE/"sensitivity_by_category_celltype.csv",index=False)
    cov=(df[df.gene=="GATA6"].groupby(["disease","cell_type"],as_index=False)
       .agg(n_sections=("sample","nunique"),min_n_cells=("n_cells","min"),
            median_n_cells=("n_cells","median"),max_n_cells=("n_cells","max"),
            n_sections_100plus=("n_cells",lambda x:int((x>=100).sum()))))
    cov.to_csv(BASE/"celltype_section_coverage.csv",index=False)
    print("Metadata verified: 5 Healthy, 9 Nonlesional, 7 Comedo, 3 Pustule",flush=True)
    print(result[(result.gene.isin(["GATA6","FABP5"])) &
          (result.cell_type.isin(["HF KC","Basal KC","Transitional KC"])) &
          (result["filter"].isin(["all","recent_min100"]))].to_string(index=False))
    print("NOTE: without confirmed donor IDs, these are descriptive section summaries, not donor-level inference.")
if __name__=="__main__":main()
