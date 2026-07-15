import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolfiles

# CONFIG
BASE_FOLDER = "data/ValidationEcoli_sansKiKd"
CSV_PATH = "/pasteur/helix/projects/Imod-moulin/Sixte/Code/Récupération_base_de_données_Chembl/Résultats/Résultats_finaux_propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Biolip2.csv"

PDB_COLUMN = "PDB_ID"
LIGAND_COLUMN = "Ligand_ID"
PROTEIN_COLUMN = "UniProt_ID"

df = pd.read_csv(CSV_PATH)

#mapping PDB → ligand ID
pdb_to_ligand = {
    str(row[PDB_COLUMN]).strip().upper(): str(row[LIGAND_COLUMN]).strip()
    for _, row in df.iterrows()
}

# MAIN

for protein_id in os.listdir(BASE_FOLDER):
    protein_dir = os.path.join(BASE_FOLDER, protein_id)

    if not os.path.isdir(protein_dir):
        continue

    ligands_dir = os.path.join(protein_dir, "ligands")
    if not os.path.exists(ligands_dir):
        continue

    ism_file = os.path.join(protein_dir, "actives_final.ism")
    open(ism_file, "w").close()

    print(f"\n🧬 Protéine {protein_id}")

    for file in os.listdir(ligands_dir):
        if not file.endswith(".pdb"):
            continue

        ligand_path = os.path.join(ligands_dir, file)

        # récupérer PDB_ID depuis le nom fichier
        pdb_id = file.replace("ligand_", "").replace(".pdb", "")

        ligand_id = pdb_to_ligand.get(pdb_id, "UNK")

        # lecture molécule
        mol = rdmolfiles.MolFromPDBFile(ligand_path, sanitize=True, removeHs=True)

        if mol is None:
            print(f"⚠️ Impossible de lire {file}")
            continue

        try:
            smi = Chem.MolToSmiles(mol)
        except:
            print(f"⚠️ SMILES échoué pour {file}")
            continue

        # écriture
        with open(ism_file, "a") as f_out:
            f_out.write(f"{smi}\t{ligand_id}\n")

        print(f"{file} ajouté")

print("\n Tous les SMILES ont été générés !")
