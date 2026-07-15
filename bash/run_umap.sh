#!/bin/bash
#SBATCH --job-name=umap_drugclip
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -o logs/umap_%x_%j.out
#SBATCH -e logs/umap_%x_%j.err

# ══════════════════════════════════════════════
# Use :
#   sbatch run_umap.sh 10
#   sbatch run_umap.sh 5 prot1 prot2 prot3 prot4 prot5
#
# Argument 1 : N (number of proteins)
# Arguments 2 optional : name of proteins
# ══════════════════════════════════════════════

BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"

source ${BASE_DIR}/bin/activate
cd ${BASE_DIR}

mkdir -p logs

N=${1:-10}
shift
TARGETS="$@"

# Case 1 : Manually described proteins
if [ -n "${TARGETS}" ]; then
    python visualize_umap_N.py \
        --n        ${N} \
        --emb_dir  ./résultats/embeddings \
        --targets  ${TARGETS} \
        --output   umap_N${N}_manuel.png \
        --seed     42

# Case 2 : N random proteins
else
    python visualize_umap_N.py \
        --n        ${N} \
        --emb_dir  ./résultats/embeddings \
        --scores_dir ./résultats \
        --targets P0A6K3 P02924
        --output   umap_N${N}_aleatoire.png \
        --seed     42
fi

