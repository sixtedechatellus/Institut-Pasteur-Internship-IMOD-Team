import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolfiles

# ===============================
# CONFIG
# ===============================
CSV_PATH = "Validationdatasets/BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv"  # ton CSV
BASE_FOLDER = "data/ValidationEcoli"                      # dossier contenant les sous-dossiers PDB_ID/
PDB_COLUMN = "PDB_ID"                                      # colonne contenant le PDB du complexe
LIGAND_COLUMN = "Ligand_ID"                                # colonne contenant l'ID du ligand

# ===============================
# MAIN
# ===============================
df = pd.read_csv(CSV_PATH)

for idx, row in df.iterrows():
    pdb_id = str(row[PDB_COLUMN]).strip().upper()
    ligand_id = str(row[LIGAND_COLUMN]).strip()
    
    # dossier du PDB
    pdb_dir = os.path.join(BASE_FOLDER, pdb_id)
    os.makedirs(pdb_dir, exist_ok=True)
    
    # fichier ligand généré précédemment
    ligand_file = os.path.join(pdb_dir, f"ligand_{pdb_id}.pdb")
    
    if not os.path.exists(ligand_file):
        print(f"⚠️ Ligand pour le complexe {pdb_id} non trouvé, skip")
        continue
    
    # Convertir le PDB du ligand en Mol
    mol = rdmolfiles.MolFromPDBFile(ligand_file, sanitize=True, removeHs=True)
    if mol is None:
        print(f"⚠️ Impossible de lire le ligand {ligand_file}, skip")
        continue
    
    # SMILES
    smi = Chem.MolToSmiles(mol)
    
    # écrire dans actives_final.ism dans le dossier du complexe
    ism_file = os.path.join(pdb_dir, "actives_final.ism")
    with open(ism_file, "a") as f_out:
        f_out.write(f"{smi}\t{ligand_id}\n")
    
    print(f"✅ Ligand {ligand_id} du complexe {pdb_id} ajouté dans {ism_file}")

print("\n✅ Tous les ligands ont été traités !")
