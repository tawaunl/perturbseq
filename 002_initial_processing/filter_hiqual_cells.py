import argparse

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc

import muon

parser = argparse.ArgumentParser()
parser.add_argument('in_path')
parser.add_argument('out_path')
parser.add_argument('cc_genes_path')

args = parser.parse_args()

#in_path = "../scratch/output_3/update_adt_umaps.h5mu"

in_path = args.in_path
out_path = args.out_path

mdata = muon.read(in_path)

keep_cells = mdata['GEX'].obs['qc_cluster'] == 'qc_clust1'
keep_cells = mdata['GEX'].obs.index[keep_cells]

# Workaround: Delete connectivities/distances, because current version
# of mudata has a bug when copying a subsetted mudata:
#
# "ValueError: Value passed for key 'connectivities' is of incorrect
# shape. Values of obsp must match dimensions ('obs', 'obs') of
# parent. Value had shape (54603, 189448) while it should have had
# (54603, 54603)."

del mdata.obsp['connectivities']
del mdata['GEX'].obsp['connectivities']
del mdata['ADT'].obsp['connectivities']
del mdata['HTO'].obsp['connectivities']

del mdata.obsp['distances']
del mdata['GEX'].obsp['distances']
del mdata['ADT'].obsp['distances']
del mdata['HTO'].obsp['distances']

mdata = mdata[keep_cells,:].copy()

# paga can crash if leiden_colors was already in uns. Also delete
# other leiden related variables to be safe
for k in ['GEX', 'ADT']:
    del mdata[k].uns['leiden']
    del mdata[k].uns['leiden_colors']
    del mdata[k].uns['leiden_sizes']

# RNA clustering

sc.pp.highly_variable_genes(mdata["GEX"])
sc.tl.pca(mdata["GEX"], use_highly_variable=True)
sc.pp.neighbors(mdata["GEX"])
sc.tl.leiden(mdata["GEX"])
sc.tl.paga(mdata["GEX"])

## Fake plot of paga to avoid error:
## ValueError: Plot PAGA first, so that adata.uns['paga']with key 'pos'.
sc.pl.paga(mdata["GEX"], show=False)

sc.tl.umap(mdata["GEX"], init_pos="paga")

# ADT clustering
sc.tl.pca(mdata['ADT'], use_highly_variable=False)
sc.pp.neighbors(mdata['ADT'])
sc.tl.leiden(mdata['ADT'])
sc.tl.paga(mdata['ADT'])

## Fake plot of paga to avoid error:
## ValueError: Plot PAGA first, so that adata.uns['paga']with key 'pos'.
sc.pl.paga(mdata["ADT"], show=False)

sc.tl.umap(mdata["ADT"], init_pos="paga")

# cell cycle

s_genes = []
g2m_genes = []
for _, row in pd.read_csv(args.cc_genes_path).iterrows():
    g = row['gene']
    g = g[0].upper() + g[1:].lower()
    if g not in mdata['GEX'].var.index:
        continue
    elif row['phase'] == "S":
        s_genes.append(g)
    elif row['phase'] == "G2M":
        g2m_genes.append(g)
    else:
        assert False

sc.tl.score_genes_cell_cycle(mdata['GEX'], s_genes, g2m_genes)

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
