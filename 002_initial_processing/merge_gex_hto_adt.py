import logging
import os
import re

import numpy as np
import pandas as pd

import anndata as ad
import pegasusio as pgio
import demuxEM
import muon

from rds2py import read_rds

import argparse

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('data_dir')
parser.add_argument('lib_csv')
parser.add_argument('dropletqc_rds')
parser.add_argument('sunrise_cellbender_gex_dir')
parser.add_argument('out_h5mu')

args = parser.parse_args()
data_dir = args.data_dir
lib_csv = args.lib_csv
dropletqc_rds = args.dropletqc_rds
out_h5mu = args.out_h5mu

ncores = int(os.environ["SLURM_CPUS_ON_NODE"])

lib_df = pd.read_csv(lib_csv)

dropletqc_res = read_rds(dropletqc_rds)
dropletqc_res = dict(zip(
    dropletqc_res['attributes']['names']['data'],
    [pd.Series(dropletqc_res['data'][i]['data'][0]['data'],
               index=dropletqc_res['data'][i]['attributes']['row.names']['data'])
     for i in range(len(dropletqc_res['data']))]
))

adata_gex_unfilt_list = []
adata_gex_cellbender_list = []
adata_hto_list = []
adata_adt_list = []

uniq_sam_ids = set(lib_df['SAMID'])
for sam_id in uniq_sam_ids:
    logging.info(sam_id)

    lib_in_samid = lib_df["SAMID"] == sam_id
    samid_lib_df = lib_df[lib_in_samid].copy().reset_index(inplace=False, drop=True)
    assert samid_lib_df.shape[0] == 5

    is_gex = samid_lib_df['Library Subtype'] == 'Expression V3.1'
    assert is_gex.sum() == 1

    is_hto = samid_lib_df['Library Subtype'].str.startswith('HTO')
    assert is_hto.sum() == 1

    is_adt = samid_lib_df['Library Subtype'] == 'ADT V3.1'
    assert is_adt.sum() == 1

    is_sgrna = samid_lib_df['Library Subtype'] == 'gRNA V3.1'
    assert is_sgrna.sum() == 2

    gex_lib, = samid_lib_df[is_gex]["Library Name"]
    adt_lib, = samid_lib_df[is_adt]["Library Name"]
    hto_lib, = samid_lib_df[is_hto]["Library Name"]

    gex_cellbender_path = os.path.join(
        args.sunrise_cellbender_gex_dir,
        f"{sam_id}/croo_output/{sam_id}_out_filtered.h5"   
    )

    gex_unfilt_path = os.path.join(data_dir, f"{gex_lib}_{sam_id}",
                                   "filtered_feature_bc_matrix.h5")

    adt_path = os.path.join(data_dir, f"{adt_lib}_{sam_id}",
                            f"{adt_lib}_{sam_id}.csv")

    hto_path = os.path.join(data_dir, f"{hto_lib}_{sam_id}",
                            f"{hto_lib}_{sam_id}.csv")

    logging.info(f"GEX unfiltered path: {gex_unfilt_path}")
    logging.info(f"GEX cellbender path: {gex_cellbender_path}")
    logging.info(f"ADT path: {adt_path}")
    logging.info(f"HTO path: {hto_path}")

    adata_gex_unfilt = pgio.read_input(gex_unfilt_path)
    adata_gex_cellbender = pgio.read_input(gex_cellbender_path)

    adata_hto = pgio.read_input(hto_path)

    # Fix mangled HTO names (strip carriage returns)
    adata_hto.var.index = [x.strip() for x in adata_hto.var.index]

    demuxEM.estimate_background_probs(adata_hto)
    demuxEM.demultiplex(adata_gex_cellbender, adata_hto, n_threads=ncores)

    adata_gex_cellbender = adata_gex_cellbender.to_anndata()
    adata_gex_unfilt = adata_gex_unfilt.to_anndata()

    adata_hto = adata_hto.to_anndata()

    adata_adt = pgio.read_input(adt_path).to_anndata()

    # Get nuclear fraction from DropletQC, and if necessary, intersect
    # it. (It has the exact same cells as the original cellranger
    # filtering, but may differ from cellbender filtering)
    nuc_frac = dropletqc_res[sam_id].copy()
    barcode_re = re.compile(r'([ACTG]+)-1')
    nuc_frac.index = [barcode_re.match(idx).group(1) for idx in nuc_frac.index]

    logging.info(f'Unfiltered GEX cells: {adata_gex_unfilt.shape[0]}')
    logging.info(f'Cellbender GEX cells: {adata_gex_cellbender.shape[0]}')

    keep = (adata_gex_unfilt.obs.index
            .intersection(adata_gex_cellbender.obs.index))

    logging.info(f"Intersection cells: {len(keep)}")

    keep = adata_hto.obs.index.intersection(keep)
    keep = adata_adt.obs.index.intersection(keep)

    logging.info(f"ADT cells: {adata_adt.shape[0]}")
    logging.info(f"HTO cells: {adata_hto.shape[0]}")
    logging.info(f"Cells after intersecting HTO and ADT: {len(keep)}")

    logging.info(f"DropletQC cells: {nuc_frac.shape[0]}")
    keep = nuc_frac.index.intersection(keep)
    logging.info(f'Cells after intersecting DropletQC: {len(keep)}')

    if len(keep) != adata_gex_cellbender.shape[0]:
        logging.warning("Intersection cells less than cellbender cells.")

    adata_gex_unfilt = adata_gex_unfilt[keep,].copy()
    adata_gex_cellbender = adata_gex_cellbender[keep,].copy()
    adata_hto = adata_hto[keep,].copy()
    adata_adt = adata_adt[keep,].copy()

    adata_gex_cellbender.obs["SAMID"] = sam_id
    adata_gex_cellbender.obs['nuclear_fraction'] = nuc_frac[adata_gex_cellbender.obs.index]

    idx = [f'{sam_id}_{i}' for i in adata_gex_cellbender.obs.index]
    adata_gex_cellbender.obs.index = list(idx)
    adata_gex_unfilt.obs.index = list(idx)
    adata_hto.obs.index = list(idx)
    adata_adt.obs.index = list(idx)

    adata_gex_cellbender_list.append(adata_gex_cellbender)
    adata_gex_unfilt_list.append(adata_gex_unfilt)
    adata_hto_list.append(adata_hto)
    adata_adt_list.append(adata_adt)

adata_gex_cellbender_concat = ad.concat(adata_gex_cellbender_list)
adata_gex_unfilt_concat = ad.concat(adata_gex_unfilt_list)
adata_hto_concat = ad.concat(adata_hto_list)
adata_adt_concat = ad.concat(adata_adt_list)

assert all(
    adata.shape[1] == adata_gex_cellbender_concat.shape[1]
    for adata in adata_gex_cellbender_list
)

assert all(
    adata.shape[1] == adata_gex_unfilt_concat.shape[1]
    for adata in adata_gex_unfilt_list
)

assert all(
    adata.shape[1] == adata_hto_concat.shape[1]
    for adata in adata_hto_list
)

assert all(
    adata.shape[1] == adata_adt_concat.shape[1]
    for adata in adata_adt_list
)

# MuData var names are global, so make sure there are no collisions
assert len(adata_gex_cellbender_concat.var.index.intersection(adata_adt_concat.var.index)) == 0
assert len(adata_gex_cellbender_concat.var.index.intersection(adata_hto_concat.var.index)) == 0
assert len(adata_adt_concat.var.index.intersection(adata_hto_concat.var.index)) == 0

mdata = muon.MuData({"GEX": adata_gex_cellbender_concat,
                     "HTO": adata_hto_concat,
                     "ADT": adata_adt_concat})

# for convenience, add the raw attribute now
for k in ['GEX', 'HTO', 'ADT']:
    mdata[k].raw = mdata[k]

# Since MuData var names are global, we make the unfiltered counts a
# layer of the cellbender anndata. But note this prevents backed
# loading of the MuData
assert all(adata_gex_unfilt_concat.obs.index == mdata['GEX'].obs.index)
assert all(adata_gex_unfilt_concat.var.index == mdata['GEX'].var.index)

mdata['GEX'].layers['uncleaned_counts'] = adata_gex_unfilt_concat.X

# add sample name

samid2sample = {}
for _, row in lib_df.iterrows():
    sam_id = row['SAMID']
    sam_name = row['Sample Name']
    if sam_id in samid2sample:
        assert samid2sample[sam_id] == sam_name
    samid2sample[sam_id] = sam_name

mdata['GEX'].obs['sample_name'] = [
    samid2sample[s] for s in mdata['GEX'].obs['SAMID']
]

# save it

mdata.write(out_h5mu)
