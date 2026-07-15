#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o logs/nettoyage_decoys_%j.out
#SBATCH -e logs/nettoyage_decoys_%j.err

# ----------------------------------------------------------
# Use :
#   sbatch script_clean_decoys.sh <FILE>
# Example :
#   sbatch script_clean_decoys.sh data/ValidationEcoli/
# ----------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch $0 <dossier_contenant_les_decoys>"
    exit 1
fi

ROOT=$1

python script_clean_decoys.py "${ROOT}"

