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

# Configurations utilisateur
BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"
DATA_DIR="$BASE_DIR/data/ValidationEcoli"   # dossier contenant les sous-dossiers PDB
DEEPCOY="$BASE_DIR/deepcoy/DeepCoy.py"
MODEL="$BASE_DIR/deepcoy/DeepCoy_DUDE_model_e09.pickle"
BATCH_SIZE=1
NUMBER_DEC=50   # nombre de décoys candidats par actif
DATASET="zinc"

# Activer l'environnement virtuel
source "$BASE_DIR/bin/activate"
echo "[✔] Environnement virtuel activé : $VIRTUAL_ENV"

for d in "$DATA_DIR"/*/ ; do
    PDB_ID=$(basename "$d")
    JSON_FILE="$d/molecules_${PDB_ID}.json"
    OUTPUT_FILE="$d/decoys_generated.smi"

    echo "[⚙] Fix v_to_keep pour DeepCoy sur $JSON_FILE"

    python -c "
import json
f = '$JSON_FILE'
data = json.load(open(f))
new_data = []
for mol in data:
    nodes = mol.get('node_features_in', [])
    if not nodes:
        print('⚠ WARNING: molécule vide ignorée')
        continue
    mol['v_to_keep'] = list(range(len(nodes)))
    new_data.append(mol)
json.dump(new_data, open(f, 'w'))
print('Fix appliqué : v_to_keep rempli pour', f, ', molécules conservées :', len(new_data))
"

    echo "[⚙] Lancement DeepCoy pour $PDB_ID"

    python "$DEEPCOY" \
        --restore "$MODEL" \
        --dataset "$DATASET" \
        --config "{\"generation\": true, \
                   \"number_of_generation_per_valid\": $NUMBER_DEC, \
                   \"batch_size\": $BATCH_SIZE, \
                   \"train_file\": \"$JSON_FILE\", \
                   \"valid_file\": \"$JSON_FILE\", \
                   \"output_name\": \"$OUTPUT_FILE\"}"

    if [ -f "$OUTPUT_FILE" ]; then
        echo "[✔] Génération terminée pour $PDB_ID -> $OUTPUT_FILE"
    else
        echo "[✖] Échec génération décoys pour $PDB_ID"
    fi
done
echo "[✔] Tous les traitements terminés."
