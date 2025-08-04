import argparse

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc

import muon

parser = argparse.ArgumentParser()
parser.add_argument('in_path')
parser.add_argument('out_path')

args = parser.parse_args()

in_path = args.in_path
out_path = args.out_path

mdata = muon.read(in_path)

# rerun ADT clustering
sc.tl.pca(mdata['ADT'], use_highly_variable=False)
sc.pp.neighbors(mdata['ADT'])
sc.tl.leiden(mdata['ADT'])
sc.tl.paga(mdata['ADT'])

## Fake plot of paga to avoid error:
## ValueError: Plot PAGA first, so that adata.uns['paga']with key 'pos'.
sc.pl.paga(mdata["ADT"], show=False)

sc.tl.umap(mdata["ADT"], init_pos="paga")

mdata_gex_adt = muon.MuData({
    "GEX": mdata["GEX"],
    "ADT": mdata["ADT"],
})

muon.pp.neighbors(mdata_gex_adt)

muon.tl.umap(mdata_gex_adt)

assert (mdata_gex_adt.obs.index == mdata['GEX'].obs.index).all()

mdata.obsm = mdata_gex_adt.obsm
mdata.obsp = mdata_gex_adt.obsp
mdata.uns = mdata_gex_adt.uns

mdata.write(out_path)
