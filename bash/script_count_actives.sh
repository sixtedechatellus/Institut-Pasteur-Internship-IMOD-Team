#!/bin/bash
#SBATCH --job-name=mean_actives
#SBATCH --output=logs/mean_actives_%j.out
#SBATCH --error=logs/mean_actives_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# ---------------------------------------------
# Use :
#   sbatch run_compute_mean_actives.sbatch <DOSSIER_RACINE>
# Example :
#   sbatch run_compute_mean_actives.sbatch data/ValidationEcoli/
# ---------------------------------------------

if [ "$#" -ne 1 ]; then
    exit 1
fi
ROOT=$1
source ~/anaconda3/etc/profile.d/conda.sh

python script_count_actives.py

