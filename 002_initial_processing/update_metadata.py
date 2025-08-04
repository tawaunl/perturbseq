import argparse
import re
import math

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc

import muon

parser = argparse.ArgumentParser()
parser.add_argument('in_path')
parser.add_argument('out_path')

args = parser.parse_args()

#in_path = "../scratch/output_3/update_adt_umaps.h5mu"

in_path = args.in_path
out_path = args.out_path

mdata = muon.read(in_path)

# Add gate

samples = set(mdata['GEX'].obs['sample_name'])

sample2gate = {}
for k in samples:
    if "High" in k:
        sample2gate[k] = "Trem2Hi"
    elif "Low" in k:
        sample2gate[k] = "Trem2Lo"
    else:
        assert "ALL" in k
        sample2gate[k] = "All"

mdata['GEX'].obs['gate'] = [sample2gate[k] for k in mdata['GEX'].obs['sample_name']]

# Add sgRNA var: gene, guide_num, guide_set
# (guide_set splits the NTC into control and dummies)

sgrna_re = re.compile(r'(.*)_(\d+)')
gene = []
guide_num = []
perturbation = []

for grna in mdata['sgRNA'].var.index:
    matched = sgrna_re.match(grna)
    curr_gene = matched.group(1)
    curr_num = int(matched.group(2))

    gene.append(curr_gene)
    guide_num.append(curr_num)

    if curr_gene != 'NTC':
        perturbation.append(curr_gene)
    else:
        if curr_num <= 100:
            perturbation.append("NTC_1_100")
        else:
            perturbation.append("NTC_{}_{}".format(
                math.floor((curr_num - 1) / 4) * 4 + 1,
                math.floor((curr_num - 1) / 4) * 4 + 4
            ))

mdata['sgRNA'].var['guide_num'] = guide_num
mdata['sgRNA'].var['guide_set'] = gene
mdata['sgRNA'].var['guide_set_with_ntc_holdouts'] = perturbation

mdata.write(out_path)
