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

#in_path = "scratch/231005_merge_gex_hto_adt_small.h5mu"
#in_path = "scratch/231005_merge_gex_hto_adt.h5mu"
in_path = args.in_path

#out_path = "scratch/231005_add_basic_stats.h5mu"
out_path = args.out_path

mdata = muon.read(in_path)

# Add QC metrics to GEX, ADT, HTO
mdata["GEX"].var['mt'] = mdata["GEX"].var_names.str.startswith('mt-')
mdata["GEX"].var['rp'] = mdata["GEX"].var_names.str.startswith('Rpl') | mdata["GEX"].var_names.str.startswith('Rps')

sc.pp.calculate_qc_metrics(
    mdata["GEX"], qc_vars=['mt', 'rp'],
    percent_top=None, inplace=True
)

sc.pp.calculate_qc_metrics(
    mdata["GEX"], qc_vars=['mt', 'rp'],
    expr_type='uncleaned_counts',
    layer='uncleaned_counts',
    percent_top=None, inplace=True
)

sc.pp.calculate_qc_metrics(mdata["ADT"], percent_top=None, inplace=True)

sc.pp.calculate_qc_metrics(mdata["HTO"], percent_top=None, inplace=True)

# Add column for n htos detected

n_htos_detected = []

for x in mdata['GEX'].obs['assignment'].astype(str):
    if x.strip() == '':
        n_htos_detected.append(0)
    else:
        n_htos_detected.append(len(x.split(',')))

mdata['GEX'].obs['n_htos_detected'] = n_htos_detected

# Add column for hto diversity

hto_diversity = np.asarray(mdata["HTO"].raw.X.todense())
n_hto = np.squeeze(hto_diversity.sum(axis=1))

#hto_diversity = (hto_diversity * (hto_diversity - 1)).sum(axis=1)
hto_diversity = (hto_diversity * hto_diversity).sum(axis=1)
hto_diveristy = np.squeeze(hto_diversity)

#hto_diversity = hto_diversity / n_hto / (n_hto-1)
hto_diversity = hto_diversity / n_hto / n_hto

hto_diversity = 1 / hto_diversity

mdata["GEX"].obs["hto_diversity"] = hto_diversity

# Normalize and log-transform

sc.pp.normalize_total(mdata["GEX"])
sc.pp.log1p(mdata["GEX"])

## We will use a separate median-polish approach (in R) to ADT
## normalization and clustering, after QC step
#sc.pp.normalize_total(mdata["ADT"])
#sc.pp.log1p(mdata["ADT"])

sc.pp.normalize_total(mdata["HTO"])
sc.pp.log1p(mdata["HTO"])

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

# HTO clustering
sc.pp.neighbors(mdata["HTO"], use_rep="X")
sc.tl.umap(mdata["HTO"])

## We will use a separate median-polish approach (in R) to ADT
## normalization and clustering, after QC step
##
## ADT clustering
#sc.tl.pca(mdata['ADT'], use_highly_variable=False)
#sc.pp.neighbors(mdata['ADT'])
#sc.tl.leiden(mdata['ADT'])
#sc.tl.paga(mdata['ADT'])
#
### Fake plot of paga to avoid error:
### ValueError: Plot PAGA first, so that adata.uns['paga']with key 'pos'.
#sc.pl.paga(mdata["ADT"], show=False)
#
#sc.tl.umap(mdata["ADT"], init_pos="paga")

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

# Save
mdata.write(out_path)
