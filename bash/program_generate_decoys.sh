#!/bin/bash
#SBATCH --job-name=decoys
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o logs/decoys_%j.out
#SBATCH -e logs/decoys_%j.err

module purge

BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"
DATA_DIR="$BASE_DIR/data/ValidationEcoli_sansKiKd"
DEEPCOY="$BASE_DIR/deepcoy/DeepCoy.py"
MODEL="$BASE_DIR/deepcoy/DeepCoy_DUDE_model_e09.pickle"
BATCH_SIZE=1
NUMBER_DEC=100
DATASET="zinc"

source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/miniconda3/etc/profile.d/conda.sh
conda activate DeepCoy-env
python -c "import tensorflow as tf; print(tf.__version__); print(tf.ConfigProto)"

cd "$BASE_DIR" || exit 1
export PYTHONPATH=$PYTHONPATH:$(pwd)

# --- Vérification des arguments DeepCoy ---
for d in "$DATA_DIR"/*/; do
    PDB_ID=$(basename "$d")
    ACTIVES_SMI="${d}actives_final.ism"
    JSON_FILE="${d}molecules_${PDB_ID}.json"
    OUTPUT_BASE="${d}decoys_generated"
    OUTPUT_FILE="${OUTPUT_BASE}.smi"
    LOG_FILE="${d}deepcoy_log_${PDB_ID}.txt"

    echo "Traitement : $PDB_ID"

    if [ ! -f "$ACTIVES_SMI" ]; then
        continue
    fi

    # STEP 1 : Preprocessing

    cd "$d" || continue

    python - <<PYEOF
import sys, os
sys.path.append("$BASE_DIR/deepcoy/data")
from prepare_data import read_file, preprocess

raw_data = read_file("actives_final.ism")
print(f"    Molécules lues : {len(raw_data)}")
preprocess(raw_data, "$DATASET", "$PDB_ID")
PYEOF

    cd "$BASE_DIR" || exit 1

    if [ ! -f "$JSON_FILE" ]; then
        find "$d" -name "*.json" 2>/dev/null
        continue
    fi

    python - <<PYEOF
import json

f = "$JSON_FILE"
data = json.load(open(f))
for i, mol in enumerate(data):
    nodes    = mol.get('node_features_in', [])
    smiles   = mol.get('smiles', 'N/A')
    v_to_keep = mol.get('v_to_keep', 'ABSENT')
PYEOF

    # STEP 3 : Fix v_to_keep

    python - <<PYEOF
import json, sys

f = "$JSON_FILE"
data  = json.load(open(f))
new_data = []
skip = 0
for mol in data:
    nodes = mol.get('node_features_in', [])
    if not nodes:
        skip += 1
        continue
    if 'v_to_keep' not in mol:
        mol['v_to_keep'] = list(range(len(nodes)))
    new_data.append(mol)

json.dump(new_data, open(f, 'w'))

if len(new_data) == 0:
    sys.exit(1)
PYEOF

    if [ $? -ne 0 ]; then
        continue
    fi

    # STEP 4 : Config et génération DeepCoy

    CONFIG_FILE="${d}config_${PDB_ID}.json"

    python - <<PYEOF
import json
config = {
    "generation": True,
    "number_of_generation_per_valid": $NUMBER_DEC,
    "batch_size": $BATCH_SIZE,
    "train_file": "$JSON_FILE",
    "valid_file": "$JSON_FILE",
    "output_name": "$OUTPUT_BASE"
}
json.dump(config, open("$CONFIG_FILE", "w"), indent=2)
print("    Config :")
print(json.dumps(config, indent=4))
PYEOF

    echo "[⚙] Lancement DeepCoy..."
    python "$DEEPCOY" \
        --restore "$MODEL" \
        --dataset "$DATASET" \
        --config "$(cat $CONFIG_FILE)" 2>&1 | tee "$LOG_FILE"


    if [ -f "$OUTPUT_FILE" ]; then
        echo "Succès : $OUTPUT_FILE ($(wc -l < $OUTPUT_FILE) lignes)"
    else
        echo "Échec génération pour $PDB_ID"
        echo "    Fichiers présents dans $d :"
        ls -lh "$d"
    fi

    rm -f "$CONFIG_FILE"
    echo
done

echo "Pipeline terminé."
