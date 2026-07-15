#!/bin/bash
#SBATCH --job-name=decoys
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -o logs/decoys_%j.out
#SBATCH -e logs/decoys_%j.err

# ============================================================
# Script SLURM : programmedecoys.sh
# Objectif : Préprocessing + génération des décoys DeepCoy
# ============================================================

module purge
module load cuda/11.8



# Répertoires et fichiers
BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"
DATA_DIR="$BASE_DIR/data/ValidationEcoli"   # dossiers PDB
DEEPCOY="$BASE_DIR/deepcoy/DeepCoy.py"
MODEL="$BASE_DIR/deepcoy/DeepCoy_DUDE_model_e09.pickle"
BATCH_SIZE=1
NUMBER_DEC=50      # nombre de décoys par actif
DATASET="zinc"

# Activer l'environnement virtuel
source "$BASE_DIR/bin/activate"
echo "[✔] Environnement virtuel activé : $VIRTUAL_ENV"

# Parcours des dossiers PDB
for d in "$DATA_DIR"/*/ ; do
    PDB_ID=$(basename "$d")
    SMILES_FILE="$d/actives_final.ism"
    
    # Vérifier que le fichier existe
    if [ ! -f "$SMILES_FILE" ]; then
        echo "[⚠] Pas d'actives_clean.smi dans $d, skip."
        continue
    fi

    JSON_FILE="$d/molecules_${PDB_ID}.json"
    OUTPUT_FILE="$d/decoys_generated.smi"

    echo "[⚙] Création JSON pour $PDB_ID avec prepare_data"

    # Appel à prepare_data / preprocess robuste
    python - <<EOF
import sys, os
sys.path.append("$BASE_DIR/deepcoy/data")
from prepare_data2 import preprocess, read_file, dataset_info

dataset = "$DATASET"
raw_data = read_file("$SMILES_FILE", reverse=True)
preprocess(raw_data, dataset, "$PDB_ID", save_dir="$d/")
EOF

    echo "[⚙] Lancement DeepCoy pour $PDB_ID"

    python "$DEEPCOY" \
        --restore "$MODEL" \
        --dataset "$DATASET" \
        --config "{\"generation\": true, \"number_of_generation_per_valid\": $NUMBER_DEC, \"batch_size\": $BATCH_SIZE}" \
        --train_file "$JSON_FILE" \
        --out_file "$OUTPUT_FILE"

done
