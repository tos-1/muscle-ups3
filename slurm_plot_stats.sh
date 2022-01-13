#!/bin/bash
#SBATCH --job-name=plot_stats
#SBATCH --output=job_stats.out
#SBATCH --error=job_stats.err
#SBATCH --mem=16000MB

source ../venv/bin/activate
srun python plot_stats.py 
deactivate

