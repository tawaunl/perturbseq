import argparse

import numpy as np
import pandas as pd

import muon
import sklearn.preprocessing
import sklearn.cluster

parser = argparse.ArgumentParser()
parser.add_argument('in_mdata')
parser.add_argument('out_mdata')
parser.add_argument('--mt_manual_cutoff', type = int, required=True)
parser.add_argument('--gene_manual_cutoff', type = int, required=True)
parser.add_argument('--hto_diversity_cutoff', type = float, required=True)
parser.add_argument('--gex_manual_cutoff', type = int, required=True)
parser.add_argument('--hto_manual_cutoff', type = int, required=True)
parser.add_argument('--adt_manual_cutoff', type = int, required=True)

args = parser.parse_args()

in_mdata = args.in_mdata
out_mdata = args.out_mdata

print("Reading mdata")
mdata = muon.read(in_mdata)
print("Finished reading mdata")

mt_fail_manual = mdata['GEX'].obs.pct_counts_mt > args.mt_manual_cutoff

gene_fail_manual = mdata['GEX'].obs.n_genes_by_counts < args.gene_manual_cutoff

hto_diversity_fail = mdata['GEX'].obs.hto_diversity > args.hto_diversity_cutoff

n_hto_fail = ~(mdata['GEX'].obs.n_htos_detected.isin((1,2,3)))

gex_fail_manual = mdata['GEX'].obs.total_counts < args.gex_manual_cutoff

hto_fail_manual = mdata['HTO'].obs.total_counts < args.hto_manual_cutoff

adt_fail_manual = mdata['ADT'].obs.total_counts < args.adt_manual_cutoff

pass_cutoffs = ~(
  mt_fail_manual |
    n_hto_fail |
    hto_diversity_fail |
    gene_fail_manual |
    gex_fail_manual |
    hto_fail_manual |
    adt_fail_manual
)

print("Running kmeans")
kmeans = sklearn.cluster.KMeans(n_clusters=2, random_state=0)

X = mdata['GEX'].obs[['nuclear_fraction','log1p_total_counts',
                      'log1p_n_genes_by_counts',
                      'pct_counts_mt', 'pct_counts_rp']]

X['n_htos'] = mdata['HTO'].obs['log1p_total_counts']
X['n_adt'] = mdata['ADT'].obs['log1p_total_counts']

X = pd.DataFrame(
    sklearn.preprocessing.scale(X.values),
    index=X.index,
    columns=X.columns
)

kmeans.fit(X.loc[pass_cutoffs])

qc_cluster_raw = kmeans.predict(X) + 1

print("Finished kmeans")

cluster_centers = pd.DataFrame(
    kmeans.cluster_centers_,
    columns=kmeans.feature_names_in_
)

# ensure that cluster 1 is the "good" one
if (cluster_centers['log1p_total_counts'][0] <
    cluster_centers['log1p_total_counts'][1]):
    qc_cluster_raw = 3 - qc_cluster_raw

qc_cluster = qc_cluster_raw.copy()

qc_cluster[(qc_cluster == 2) & (mdata['GEX'].obs['nuclear_fraction'] < .08)] = 3

qc_cluster = np.array(
    [f'qc_clust{i}' for i in qc_cluster], dtype='object'
)

qc_cluster[~pass_cutoffs] = 'outlier'
qc_cluster[pass_cutoffs & (mdata['GEX'].obs['n_htos_detected'] >= 2)] = 'hto_doublet'
qc_cluster[pass_cutoffs & (mdata['GEX'].obs['n_htos_detected'] == 0)] = 'hto_unassigned'

mdata['GEX'].obs['is_outlier'] = ~pass_cutoffs
mdata['GEX'].obs['hto_multiplet'] = mdata['GEX'].obs['n_htos_detected'] >= 2
mdata['GEX'].obs['pass_gex_adt_hto'] = pass_cutoffs & (mdata['GEX'].obs['n_htos_detected'] == 1)
mdata['GEX'].obs['qc_cluster'] = qc_cluster
mdata['GEX'].obs['qc_cluster_raw'] = qc_cluster_raw

print("Saving mdata")
mdata.write(out_mdata)
print("Finished")
