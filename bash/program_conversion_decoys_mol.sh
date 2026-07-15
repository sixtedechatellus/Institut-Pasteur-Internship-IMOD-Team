#!/bin/bash
#SBATCH --job-name=decoys
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH -o logs/decoys_%j.out
#SBATCH -e logs/decoys_%j.err

module purge
module load cuda/11.8

# ---------------------------
#  Activation de l’environnement Python
# ---------------------------
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate
PYTHON=$(which python)

module load openbabel/3.1.1

SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/Ecoli_UMAP"

for d in "$SRC"/*/; do
    ID=$(basename "$d")
    SMIFILE="${d}decoys_final.ism"
    MOL2_OUT="${d}decoys_final.mol2"

    if [ ! -f "$SMIFILE" ]; then
        echo "⚠️  Aucun fichier actives_clean.ism trouvé pour $ID"
        continue
    fi
    # Conversion SMILES vers MOL2
    obabel -ismi "$SMIFILE" -omol2 -O "$MOL2_OUT" 2>/dev/null

    if [ -s "$MOL2_OUT" ]; then
        echo "   ✔ $MOL2_OUT créé"
    else
        echo "   ⚠️ Échec de conversion pour $SMIFILE"
    fi
done
echo "Conversion terminée."
