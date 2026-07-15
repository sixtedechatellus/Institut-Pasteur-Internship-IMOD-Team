SRC="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip"



for d in "$SRC"/*/; do
    ID=$(basename "$d")
    PDB_FILE="${d}${ID}.pdb"
    LIGAND_PDB="${d}ligand.pdb"
    PROT_PDB="${d}receptor.pdb"
    MOL2_OUT="${d}crystal_ligand.mol2"

    echo "🔹 Traitement ${ID}"

    # Vérifie que le PDB existe
    if [ ! -f "$PDB_FILE" ]; then
        echo "   ⚠️  Fichier PDB introuvable : $PDB_FILE"
        continue
    fi

    # 1️⃣ Sépare les lignes ATOM / HETATM
    grep "^ATOM"  "$PDB_FILE" > "$PROT_PDB"
    grep "^HETATM" "$PDB_FILE" > "$LIGAND_PDB"

    # 2️⃣ Conversion ligand.pdb → crystal_ligand.mol2 via RDKit
    "$PYTHON" - <<PYEOF
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import sys

ligand_pdb = "$LIGAND_PDB"
mol2_out = "$MOL2_OUT"

try:
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    if mol is None:
        print(f"   ⚠️  Ligand vide pour {ligand_pdb}")
        sys.exit(0)
    rdmolfiles.MolToMol2File(mol, mol2_out)
    print(f"   ✔️  crystal_ligand.mol2 créé pour {ligand_pdb}")
except Exception as e:
    print(f"   ❌  Erreur conversion {ligand_pdb}: {e}")
PYEOF

done

echo "[✔] Extraction terminée."
