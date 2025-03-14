#!/bin/bash
#SBATCH --job-name=MOFA            # Job name
#SBATCH --output=output_%A_%a.out       # Output file
#SBATCH --error=error_%A_%a.err         # Error file
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks per node
#SBATCH --mem=95GB                       # Memory requested per node
#SBATCH --array=0-15                     # Array job indices (modify this)

# List of factors to iterate over (modify this to your specific needs)
factors=( $(seq 6 2 20) )
factors+=( $(seq 25 5 50) )

# Load R module if necessary
# module load R

# Get current factor based on SLURM array task ID
current_factor=${factors[$SLURM_ARRAY_TASK_ID]}

# Execute your R script with the current factor as a command line argument
Rscript /gstore/scratch/u/lucast3/ngs5388_collaborative/MOFA/MOFAIter.R "$current_factor"
