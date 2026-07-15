#!/bin/bash
#SBATCH --job-name=decoys
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH -o logs/decoys_%j.out
#SBATCH -e logs/decoys_%j.err

module purge
module load cuda/11.8

# ---------------------------
#  Activation de l’environnement Python
# ---------------------------
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate
PYTHON=$(which python)
echo "[✔] Environnement virtuel activé : $VIRTUAL_ENV"
echo "[✔] Exécutable Python : $PYTHON"

module load openbabel/3.1.1

SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip"

echo "🔹 Lancement du renommage dans $SRC"
echo "----------------------------------------------------------"

SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip"

for d in "$SRC"/*/; do
    if [ -f "${d}actives_not_clean.ism" ]; then
        cp "${d}actives_not_clean.ism" "${d}actives_final.ism"
        echo "✅  Copie créée : ${d}actives_final.ism"
    fi
done

echo "[✔]  Toutes les copies ont été créées."
