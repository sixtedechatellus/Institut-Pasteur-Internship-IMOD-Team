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

SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/Ecoli_UMAP"
module load openbabel/3.1.1 

for d in "$SRC"/*/; do
    ID=$(basename "$d")
    PDB_FILE="${d}${ID}.pdb"
    LIGAND_PDB="${d}ligand_${ID}.pdb"
    PROT_PDB="${d}receptor.pdb"
    MOL2_OUT="${d}crystal_ligand.mol2"

    echo "🔹 Traitement ${ID}"

    if [ ! -f "$PDB_FILE" ]; then
        echo "   ⚠️  Fichier PDB introuvable : $PDB_FILE"
        continue
    fi

    # 1️⃣ Sépare les lignes ATOM / HETATM
    grep "^ATOM"  "$PDB_FILE" > "$PROT_PDB"
    grep "^HETATM" "$PDB_FILE" > "$LIGAND_PDB"

# ------------------------------------------------------------
# 2️⃣ Conversion ligand.pdb → crystal_ligand.mol2 via OPEN BABEL
# ------------------------------------------------------------

obabel "$LIGAND_PDB" -O "$MOL2_OUT" -p 7.4 2>/dev/null

if [ -s "$MOL2_OUT" ]; then
    echo "   ✔ $MOL2_OUT créé"
else
    echo "   ⚠ Conversion échouée pour $LIGAND_PDB"
fi
done

echo "[✔] Extraction terminée."
