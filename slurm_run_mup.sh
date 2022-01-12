#!/bin/bash
#SBATCH --job-name=muscleups
#SBATCH --output=job_mup.out
#SBATCH --error=job_mup.err
#SBATCH --mem=16000MB
#SBATCH -n 1
#SBATCH -p squire8
# -w node8

source ../venv/bin/activate
srun python run.py -paramfile='params.py' 
deactivate

