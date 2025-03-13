#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -J 020_run_glmgampoi_singlecell_genelvl
#SBATCH --qos=medium
#SBATCH --array=1-20
#SBATCH -c 7
#SBATCH --mem=100G
#SBATCH -o log/%x_%a.out
#SBATCH -e log/%x_%a.err
#SBATCH --open-mode truncate
#SBATCH --mail-type=BEGIN,END,FAIL # Type of email notification- BEGIN,END,FAIL,ALL

# NOTE CAUTION CAUTION CAUTION!
# Make sure NJOBS matches the sbatch argument
NJOBS=20

# Following sets NJOBS automatically, but will be incorrect when
# resubmitting individual jobs!
#NJOBS="$SLURM_ARRAY_TASK_MAX"

ml R/prd

source /gstore/scratch/u/lucast3/ngs5388_collaborative/030_rna_diff_expr/paths.sh

OUTDIR="$ANALYSIS_SCRATCH"/"020_run_glmgampoi_singlecell_genelvl"


Rscript /gstore/scratch/u/lucast3/ngs5388_collaborative/004_DE_and_TopicModeling/run_glmgampoi_singlecell.R \
        --array_jobs "$NJOBS" \
        --array_idx "$SLURM_ARRAY_TASK_ID" \
        --n_subproc 7 \
        --no_shrink_overdisp \
        --ridge_prior_sigma 1 \
        "$ANALYSIS_SCRATCH"/"converted_sce.Rds" \
        "$OUTDIR"/"$SLURM_ARRAY_TASK_ID".csv
