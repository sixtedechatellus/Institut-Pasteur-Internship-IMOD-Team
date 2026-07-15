#!/bin/bash
#SBATCH --job-name=extract_embeddings
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH -o logs/embed_%x_%j.out
#SBATCH -e logs/embed_%x_%j.err

# ───────────────────────────────
# Usage : sbatch extract_embeddings.sh <nom_dataset>
# ───────────────────────────────

DATASET=$1
BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"
FOLD=0   # un seul fold suffit pour la visu

module load cuda/11.8
source ${BASE_DIR}/bin/activate

cd ${BASE_DIR}
export PYTHONPATH=$PYTHONPATH:$(pwd)

weight_path="./data/model_weights/6_folds/fold_${FOLD}.pt"
data_path="${BASE_DIR}/data/${DATASET}/"

echo "[🚀] Extraction des embeddings - $(date)"

# ── Étape 1 : inférence normale (sauvegarde automatique les .npy) ──
CUDA_VISIBLE_DEVICES="0" python ./unimol/test.py --user-dir ./unimol \
       $data_path \
       --results-path ./résultats/${DATASET} \
       --num-workers 8 --ddp-backend=c10d --batch-size 64 \
       --task drugclip --loss in_batch_softmax --arch drugclip \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256 --seed 1 \
       --use-folds False \
       --path $weight_path \
       --log-interval 100 --log-format simple \
       --max-pocket-atoms 511 \
       --test-task BIOLIP

echo "[✔] Embeddings sauvegardés dans ./résultats/embeddings/"
echo "[📊] Lancement UMAP - $(date)"

# ── Étape 2 : visualisation UMAP (CPU, pas besoin de GPU) ──
python visualize_umap.py \
    --emb_dir ./résultats/embeddings \
    --scores_dir ./résultats \
    --targets P0A6K3 P02924 \
    --output  ./résultats/umap_${DATASET}.png

echo "[✔] Terminé - $(date)"
