"""Inspect FULL GEO GSE301280 annotated h5ad, with verified cell annotations.
This is descriptive reanalysis ONLY; does not make pseudo-independent cell-level claims.
The public file is ~758 MB; no clinical or intervention data are present.
"""
import argparse
import json
import re
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

GENES=["GATA6","FABP5","CRABP2","KRT5","KRT79","KRT10","PPARG","AWAT2","RARA","RXRA","IL1B","FASN","SCD","MKI67"]
SOURCE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301280/suppl/GSE301280_anndata_global_6_29_25.h5ad"

def text_summary(s):
    x=s.astype('string').fillna("<missing>")
    counter=x.value_counts(dropna=False)
    return {"nonmissing":int((x!="<missing>").sum()),"n_unique":int(x.nunique()),
            "top_values":{str(k):int(v) for k,v in counter.head(35).items()}}

def select_field(columns, patterns):
    for pattern in patterns:
        matches=[k for k in columns if re.search(pattern,str(k),flags=re.I)]
        if matches:return matches[0]
    return None

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--h5ad",required=True)
    pa.add_argument("--out",required=True)
    pa.add_argument("--chunk-size",type=int,default=3000)
    a=pa.parse_args()
    dest=Path(a.out);dest.mkdir(parents=True,exist_ok=True)
    print("READING",a.h5ad,flush=True)
    adata=ad.read_h5ad(a.h5ad,backed="r")
    obs=adata.obs.copy()
    columns=list(obs.columns)
    var=adata.var.copy()
    checks={}
    checks["shape"]=[int(adata.n_obs),int(adata.n_vars)]
    checks["obs_columns"]=[str(x) for x in columns]
    checks["var_columns"]=[str(x) for x in var.columns]
    checks["obsm_keys"]=list(adata.obsm.keys())
    checks["layers"]=list(adata.layers.keys())
    checks["uns_keys"]=list(adata.uns.keys())
    checks["source"]=SOURCE
    # Descriptive metadata cardinality and frequency; do NOT infer donor labels from specimen prefixes.
    checks["column_summary"]={str(k):text_summary(obs[k]) for k in columns if obs[k].nunique(dropna=False)<=1000}
    obskeys=[str(k) for k in columns]
    ct=select_field(obskeys,[r"^bulk_labels$",r"cell.?type",r"cell.?label",r"annotation"])
    sample=select_field(obskeys,[r"^sample$",r"sample.?id",r"^sample_name$",r"section",r"specimen",r"^orig.ident$"])
    state=select_field(obskeys,[r"^disease_state$",r"^disease$",r"^condition$",r"^disease.*status",r"^lesion",r"^group$"])
    donor=select_field(obskeys,[r"^donor(_id)?$",r"^patient(_id)?$",r"^subject(_id)?$",r"individual",r"donor",r"patient",r"subject"])
    keys={"cell_type":ct,"sample":sample,"disease":state,"donor":donor}
    checks["selected_metadata_keys"]=keys
    (dest/"h5ad_schema.json").write_text(json.dumps(checks,indent=2,default=str))
    print("FIELDS",keys,"SHAPE",checks["shape"],flush=True)
    if ct is None:
        (dest/"BLOCKER.txt").write_text("Cell type annotation field absent. Do not interpret results.\n")
        raise RuntimeError("bulk_labels and fallback cell annotations missing")
    if sample is None:
        (dest/"BLOCKER.txt").write_text("No verified sample/section column. Gene-level analyses require mapping.\n")
        raise RuntimeError("Sample/section field missing")
    # Gene name identities verified from var before analysis, not assumed.
    feature_names=pd.Index(adata.var_names.astype(str))
    idx={g:int(np.where(feature_names==g)[0][0]) for g in GENES if np.sum(feature_names==g)==1}
    missing=[g for g in GENES if g not in idx]
    checks["present_genes"]=list(idx)
    checks["missing_genes"]=missing
    if not idx:raise RuntimeError("No target genes in annotated h5ad")
    # Keep donor blank if it is absent: sample prefixes are NOT independently confirmed IDs.
    fields={"sample":sample,"cell_type":ct}
    if state:fields["disease"]=state
    if donor:fields["donor"]=donor
    df=pd.DataFrame({new:obs[col].astype("string").fillna("<missing>").astype(str).to_numpy() for new,col in fields.items()})
    df["n_cells"]=1
    grouped=df.groupby(list(fields),dropna=False,sort=False).size().reset_index(name="n_cells")
    grouped.to_csv(dest/"annotated_cell_composition.csv",index=False)
    # Summarize cell-type labels, sample counts and disease categories, without retaining patient data.
    print("LABELS",grouped.groupby("cell_type")["n_cells"].sum().sort_values(ascending=False).head(30).to_dict(),flush=True)
    # Iterate in backed chunks so this fits modest runner memory.
    use_names=list(idx)
    indices=np.array([idx[g] for g in use_names],dtype=int)
    ng= len(use_names)
    labels=list(fields)
    acc={}
    n=adata.n_obs
    for st in range(0,n,a.chunk_size):
        en=min(st+a.chunk_size,n)
        mat=adata.X[st:en,:]
        if sparse.issparse(mat): mat=mat[:,indices].toarray()
        else: mat=np.asarray(mat)[:,indices]
        mat=np.asarray(mat,dtype=np.float64)
        sub=df.iloc[st:en]
        gidx=sub.groupby(labels,sort=False).indices
        for label,ii in gidx.items():
            if not isinstance(label,tuple): label=(label,)
            vals=mat[ii]
            if label not in acc:acc[label]=[0,np.zeros(ng),np.zeros(ng),np.zeros(ng)]
            a0=acc[label];a0[0]+=len(ii)
            a0[1]+=vals.sum(axis=0)
            a0[2]+=(vals>0).sum(axis=0)
            a0[3]+=(vals*vals).sum(axis=0)
        if st==0 or en==n or en%(a.chunk_size*25)==0:
            print("PROGRESS",en,"/",n,flush=True)
    rows=[]
    for key,(count,total,detect,sumsq) in acc.items():
        for i,gene in enumerate(use_names):
            rows.append({**dict(zip(labels,key)),"gene":gene,"n_cells":count,
                         "mean_anndata_X":total[i]/count,
                         "fraction_nonzero":detect[i]/count,
                         "variance_anndata_X":max(0.,sumsq[i]/count-(total[i]/count)**2)})
    result=pd.DataFrame(rows)
    result.to_csv(dest/"annotated_sample_celltype_gene_expression.csv",index=False)
    # Section-level comparators are used only when each section's disease category is present.
    if state:
        z=result.groupby(["disease","cell_type","gene"],as_index=False).agg(
            n_annotated_sample_celltype_groups=("sample","nunique"),
            median_sample_group_mean=("mean_anndata_X","median"),
            mean_sample_group_mean=("mean_anndata_X","mean"),
            median_fraction_nonzero=("fraction_nonzero","median"))
        z.to_csv(dest/"annotated_category_celltype_gene_summary.csv",index=False)
        print("TARGET SUMMARY",z[z.gene.isin(["GATA6","FABP5","CRABP2"])].head(80).to_string(index=False),flush=True)
    # Exact within-donor comparisons cannot be justified without a validated donor column.
    (dest/"analysis_limitations.json").write_text(json.dumps({
        "donor_metadata_verified":donor is not None,
        "donor_mapping":donor or "UNAVAILABLE: no independently verified donor key; do not infer from sample title",
        "sample_key":sample,"celltype_key":ct,"disease_key":state,
        "note":"Descriptive means of deposited annotated matrix X. X could be normalized/log transformed; review dataset provenance. Samples from same donor may be correlated. No donor-level p-values or causal conclusions.",
        "n_cells":n,"n_features":adata.n_vars,
        "missing_genes":missing},indent=2))
    print("SUCCESS: annotated cell-type summary computed",flush=True)
    adata.file.close()

if __name__=="__main__":main()
