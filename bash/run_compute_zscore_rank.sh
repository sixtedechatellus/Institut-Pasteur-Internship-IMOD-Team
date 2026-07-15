#!/bin/bash
#SBATCH --job-name=zscore_rank
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH -o logs/zscore_rank_%j.out
#SBATCH -e logs/zscore_rank_%j.err

module purge
module load cuda/11.8

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/


# ---------------------------------------------------------
# Use :
#   sbatch run_compute_zscore_rank.sbatch <DOSSIER_CSV>
# Example :
#   sbatch run_compute_zscore_rank.sbatch data/ValidationEcoli/csv_scores
# ---------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch $0 <dossier_contenant_les_CSV>"
    exit 1
fi

CSV_DIR=$1

python compute_zscore_rank.py "${CSV_DIR}"


