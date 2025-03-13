
0. Edit `config.yaml`, and set `scratch_dir` to a scratch directory
   where you have write access.

1. First, install and activate the conda env. Note that this pipeline
   uses R within conda, rather than from CEDAR, to better facilitate
   R-python interop via reticulate.

```{sh}
conda env create -f environment.yaml
conda activate ngs5388_snakemake
```

2. Within the conda R, manually install DropletQC:

```{R}
devtools::install_github("powellgenomicslab/DropletQC")
```

3. Run `snakemake`.

4. After successfully completing the pipeline, proceed onto the folder
   [003_initial_plotting_dsdb_upload](../003_initial_plotting_dsdb_upload). (Note
   the conda environment is no longer needed for this -- the standard
   CEDAR R is used for the downstream plotting and analyses after the
   initial Snakemake pipeline).

