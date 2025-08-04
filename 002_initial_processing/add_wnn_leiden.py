import argparse
import re
import math

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc

import muon

parser = argparse.ArgumentParser()
parser.add_argument('leiden_res', type=float)
parser.add_argument('leiden_res2', type=float)
parser.add_argument('in_path')
parser.add_argument('out_path')

args = parser.parse_args()

#in_path = "../scratch/output_3/filter_hiqual_cells.h5mu"

in_path = args.in_path
out_path = args.out_path
leiden_res = args.leiden_res
leiden_res2 = args.leiden_res2

mdata = muon.read(in_path)

adata_leiden = ad.AnnData(
    obs=mdata['GEX'].obs,
    uns=mdata.uns,
    obsp=mdata.obsp
)

sc.tl.leiden(adata_leiden, resolution=leiden_res)

sc.tl.leiden(adata_leiden, resolution=leiden_res2,
             key_added='leiden2')

assert all(mdata.obs.index == mdata['GEX'].obs.index)
mdata.obs['leiden'] = adata_leiden.obs['leiden']
mdata.obs['leiden2'] = adata_leiden.obs['leiden2']

mdata.write(out_path)
