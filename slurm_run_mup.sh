#!/bin/bash
#SBATCH --job-name=muscleups
#SBATCH --output=job_mup.out
#SBATCH --error=job_mup.err
# -N 1
#SBATCH -n 1
# -w node2

srun python run.py -paramfile='params.py' 

