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

SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/Ecoli_UMAP"
module load openbabel/3.1.1 

for d in "$SRC"/*/; do
    ID=$(basename "$d")
    PROT_PDB="${d}receptor.pdb"
    MOL2_OUT="${d}crystal_ligand.mol2"
    LIG_DIR="${d}ligands"

    echo "🔹 Traitement ${ID}"

    # Vérifie dossier ligands
    if [ ! -d "$LIG_DIR" ]; then
        echo "   ⚠️  Pas de dossier ligands"
        continue
    fi

    # Prend le premier ligand trouvé
    FIRST_LIGAND=$(ls "$LIG_DIR"/*.pdb 2>/dev/null | head -n 1)

    if [ -z "$FIRST_LIGAND" ]; then
        echo "   ⚠️  Aucun ligand trouvé"
        continue
    fi

    echo "   ➜ Ligand choisi : $(basename "$FIRST_LIGAND")"

    # Conversion en MOL2
    obabel "$FIRST_LIGAND" -O "$MOL2_OUT" -p 7.4 2>/dev/null

    if [ -s "$MOL2_OUT" ]; then
        echo "   ✔ $MOL2_OUT créé"
    else
        echo "   ⚠ Conversion échouée pour $FIRST_LIGAND"
    fi
done

echo "Extraction terminée."
