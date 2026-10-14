"""GSE301280: independent, descriptive Xenium matrix audit.
This script analyzes PUBLIC GEO data, but cannot infer donor identity or causal effects.
Use Python 3.11: pip install requests pandas numpy scipy h5py
"""
import argparse
import re
from pathlib import Path
from urllib.parse import unquote, urljoin
from collections import Counter
import numpy as np
import pandas as pd
import requests
import h5py
from scipy.sparse import csc_matrix
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

NAMES = """
GSM9079610 11027_A2_Healthy
GSM9079611 11027_B1_Nonlesional
GSM9079612 11027_B2_Nonlesional
GSM9079613 11027_C1_Healthy
GSM9079614 11027_C2_Healthy
GSM9079615 11027_D2_Nonlesional
GSM9079616 10305_A1_Nonlesional
GSM9079617 10305_B1_Comedo
GSM9079618 10305_C1_Comedo
GSM9079619 10305_D1_Comedo
GSM9079620 11629-JD-1_B3_Pustule
GSM9079621 11629-JD-1_A1_Nonlesional
GSM9079622 11629-JD-1_A2_Nonlesional
GSM9079623 11629-JD-1_C2_Pustule
GSM9079624 11629-JD-1_B1_Nonlesional
GSM9079625 11629-JD-2_C1_Comedone
GSM9079626 11629-JD-2_C2_Comedone
GSM9079627 11629-JD-2_D2_Comedone
GSM9079628 11629-JD-2_B2_Nonlesional
GSM9079629 11629-JD-1_C1_Comedone
GSM9079630 11629-JD-2_A2_Pustule
GSM9079631 11629-JD-2_A1_Nonlesional
GSM9079632 11629-JD-1_C3_Nonlesional
GSM9079633 11629-JD-2_D1_Healthy
""".strip()
SAMPLES = dict(line.split() for line in NAMES.splitlines())
GENES = ["GATA6","FABP5","CRABP2","KRT5","KRT79","KRT10","JUN","JUNB",
         "FOS","FOSB","PPARG","AWAT2","RARA","RARB","RXRA","RXRB",
         "MKI67","IL1B","SCD","FASN"]
CATS=["Healthy","Nonlesional","Comedo","Pustule"]

def category(label):
    for c in CATS:
        if c.lower() in label.lower(): return c
    if "comedone" in label.lower(): return "Comedo"
    raise ValueError("Unknown tissue label: "+label)

def geo_h5_url(gsm, sess):
    base=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
    r=sess.get(base,timeout=45);r.raise_for_status()
    paths=re.findall(r'href\\s*=\\s*["\\\']([^"\\\']+\\.h5(?:\\?[^"\\\']*)?)["\\\']',r.text,re.I)
    candidates=[p.split("?")[0] for p in paths if "cell_feature_matrix.h5" in p]
    if len(candidates)!=1: raise RuntimeError(f"{gsm}: expected exactly 1 H5 listing; found {candidates!r} at {base}")
    return urljoin(base,unquote(candidates[0]))

def get_file(sess, url, dest):
    if dest.exists() and dest.stat().st_size>1000:return
    tmp=dest.with_suffix(".partial")
    with sess.get(url,timeout=(25,120),stream=True) as r:
        r.raise_for_status()
        with tmp.open("wb") as f:
            for chunk in r.iter_content(262144):
                if chunk:f.write(chunk)
    if tmp.stat().st_size<1000:raise RuntimeError("H5 unexpectedly small: "+url)
    tmp.replace(dest)

def read_matrix(path):
    with h5py.File(path,"r") as h:
        m=h["matrix"]
        names=[v.decode() if isinstance(v,bytes) else str(v) for v in m["features"]["name"][()]]
        shape=tuple(map(int,m["shape"][()]))
        assert len(names)==shape[0], "Matrix orientation changed: "+str(shape)
        mat=csc_matrix((m["data"][()],m["indices"][()],m["indptr"][()]),shape=shape)
        totals=np.asarray(mat.sum(axis=0)).ravel().astype(float)
        rows=[]
        for gene in GENES:
            idx=[j for j,g in enumerate(names) if g==gene]
            if not idx:
                rows.append(dict(gene=gene,present_in_panel=False,mean_log1p_cptt=None,fraction_detected=None,mean_raw_umi=None))
                continue
            if len(idx)!=1:raise RuntimeError("Duplicate gene: "+gene)
            x=mat.getrow(idx[0]).toarray().ravel().astype(float)
            rows.append(dict(gene=gene,present_in_panel=True,
                  mean_log1p_cptt=float(np.mean(np.log1p(10000*x/np.maximum(totals,1)))),
                  fraction_detected=float(np.mean(x>0)),
                  mean_raw_umi=float(np.mean(x))))
        return shape[0],shape[1],rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="research/2026-10-14/gse301280_results")
    p.add_argument("--download",action="store_true")
    args=p.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    meta=pd.DataFrame([dict(gsm=g,title=t,category=category(t)) for g,t in SAMPLES.items()])
    meta.to_csv(out/"sample_geo_title_audit.csv",index=False)
    paper=dict(Healthy=5,Nonlesional=9,Comedo=7,Pustule=3)
    c=Counter(meta.category)
    count=pd.DataFrame([dict(category=k,geo_title_count=c[k],paper_count=paper[k],
                     geo_minus_paper=c[k]-paper[k]) for k in CATS])
    count.to_csv(out/"geo_vs_paper_category_counts.csv",index=False)
    print(count.to_string(index=False),flush=True)
    if not args.download:return
    sess=requests.Session()
    sess.headers.update({"User-Agent":"FAR2 public-acne-dataset-independent-audit/1.0"})
    sess.mount("https://",HTTPAdapter(max_retries=Retry(total=3,backoff_factor=1.5,status_forcelist=[429,500,502,503,504])))
    records=[];errors=[]
    for gsm,title in SAMPLES.items():
        try:
            url=geo_h5_url(gsm,sess)
            cache=out/"cache";cache.mkdir(exist_ok=True)
            file=cache/(gsm+"_cell_feature_matrix.h5")
            get_file(sess,url,file)
            n_genes,n_cells,stats=read_matrix(file)
            for row in stats:records.append(dict(gsm=gsm,geo_title=title,category=category(title),
                              n_genes=n_genes,n_cells=n_cells,**row))
            print(f"READ {gsm}: {n_cells} cells and {n_genes} panel features",flush=True)
        except Exception as e:
            errors.append(f"{gsm} {type(e).__name__}: {e}")
            print("ERROR",errors[-1],flush=True)
    (out/"errors.txt").write_text("\n".join(errors))
    if records:
        df=pd.DataFrame(records)
        df.to_csv(out/"sample_gene_expression.csv",index=False)
        summary=(df[df.present_in_panel]
                 .groupby(["category","gene"],as_index=False)
                 .agg(n_sections=("gsm","nunique"),median_section_expression=("mean_log1p_cptt","median"),
                      mean_section_expression=("mean_log1p_cptt","mean"),
                      median_detection_fraction=("fraction_detected","median")))
        summary.to_csv(out/"category_gene_expression.csv",index=False)
        print("PANEL PRESENCE AND SECTION-LEVEL MEDIANS; NOT INDEPENDENT DONOR TESTS",flush=True)
        for gene in ("GATA6","FABP5","CRABP2"):
            print(gene,summary[summary.gene==gene].to_dict("records"),flush=True)
    print(f"Completed {len(set(r['gsm'] for r in records))}/24 matrices; {len(errors)} errors",flush=True)
    if errors:raise RuntimeError("Dataset incomplete; no final interpretation allowed.")

if __name__=="__main__":main()
