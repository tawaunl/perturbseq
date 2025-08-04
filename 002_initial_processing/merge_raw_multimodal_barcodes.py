import logging
import os
import re

import numpy as np
import pandas as pd

import anndata as ad
import pegasusio as pgio
import demuxEM
import muon

#from rds2py import read_rds

import argparse

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('lib_parent_dir')
parser.add_argument('lib_csv')
parser.add_argument('out_path')
args = parser.parse_args()

#lib_parent_dir = "/gstore/data/ctgbioinfo/tgi_data_concierge/ngs5388/output"
#lib_df = pd.read_csv("data/Pilot 9 Trem2 FE perturb - Sample_mapping_for_MThitslims.csv")
#
#out_path = "scratch/231009_merge_raw_multimodal_barcodes.txt"

lib_parent_dir = args.lib_parent_dir
lib_df = pd.read_csv(args.lib_csv)
out_path = args.out_path

uniq_sam_ids = set(lib_df['SAMID'])

with open(out_path, "w") as f:
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
        hto_lib, = samid_lib_df[is_hto]["Library Name"]
        adt_lib, = samid_lib_df[is_adt]["Library Name"]

        adata_gex = pgio.read_input(os.path.join(lib_parent_dir, f"{gex_lib}_{sam_id}",
                                                "raw_feature_bc_matrix.h5"))
        adata_hto = pgio.read_input(os.path.join(lib_parent_dir, f"{hto_lib}_{sam_id}",
                                                f"{hto_lib}_{sam_id}.csv"))

        # Fix mangled HTO names (strip carriage returns)
        adata_hto.var.index = [x.strip() for x in adata_hto.var.index]

        adata_gex = adata_gex.to_anndata()
        adata_hto = adata_hto.to_anndata()

        adata_adt = pgio.read_input(os.path.join(lib_parent_dir, f"{adt_lib}_{sam_id}",
                                                f"{adt_lib}_{sam_id}.csv")).to_anndata()

        keep = adata_gex.obs.index.intersection(adata_hto.obs.index)
        keep = adata_adt.obs.index.intersection(keep)

        logging.info(f"ADT cells: {adata_adt.shape[0]}")
        logging.info(f"HTO cells: {adata_hto.shape[0]}")
        logging.info(f"GEX cells: {adata_gex.shape[0]}")
        logging.info(f"Intersection: {len(keep)}")

        for i in keep:
            print(f'{sam_id}_{i}', file = f)
