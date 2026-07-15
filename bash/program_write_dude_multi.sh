#!/bin/bash
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

MOL_DATA_ROOT="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/Ecoli_UMAP"
SCRIPT_PATH="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/write_dude_multi.py"

echo "[📂] Dossier racine : $MOL_DATA_ROOT"
echo "[🚀] Début du traitement de toutes les cibles..."

# ==============================================
# Boucle sur chaque sous-dossier
# ==============================================
for d in "$MOL_DATA_ROOT"/*/; do
    ID=$(basename "$d")
    echo
    echo "➡️  Traitement du dossier : $ID"
    
    # Chaque sous-dossier est traité indépendamment
    python "$SCRIPT_PATH" --mol_data_path "$d" >"${d}/run.out" 2>"${d}/run.err"
    
    STATUS=$?
    if [ $STATUS -ne 0 ]; then
        echo "⚠️  Erreur dans $ID (code=$STATUS), on continue avec la suite."
    else
        echo "✅  Terminé pour $ID"
    fi
done

echo
echo "[✔]  Tous les dossiers ont été parcourus."
