import os
import pandas as pd
import requests
from tqdm import tqdm
from Bio.PDB import PDBParser, PDBIO, Select

# ===============================
# CONFIG
# ===============================
CSV_PATH = "/pasteur/helix/projects/Imod-moulin/Sixte/Code/Récupération_base_de_données_Chembl/Résultats/Résultats_finaux_propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Biolip2.csv"
OUTPUT_DIR = "data/ValidationEcoli_sansKiKd"

PDB_COLUMN = "PDB_ID"
PROTEIN_COLUMN = "UniProt_ID"

# ===============================
# SELECTORS
# ===============================
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "

class LigandSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0].startswith("H") and residue.resname not in ["HOH", "WAT"]:
            return True
        return False

# ===============================
# DOWNLOAD
# ===============================
def download_pdb(pdb_id, dest_folder):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    dest_file = os.path.join(dest_folder, f"{pdb_id}.pdb")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and "ATOM" in r.text:
            with open(dest_file, "w") as f:
                f.write(r.text)
            return dest_file
    except:
        pass
    return None

# ===============================
# EXTRACTION
# ===============================
def save_receptor(structure, out_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ProteinSelect())

def save_ligand(structure, out_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, LigandSelect())

# ===============================
# MAIN
# ===============================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)

    grouped = df.groupby(PROTEIN_COLUMN)

    print(f"🔍 {len(grouped)} protéines uniques")

    for protein_id, group in tqdm(grouped, desc="Protéines"):

        protein_dir = os.path.join(OUTPUT_DIR, str(protein_id))
        ligands_dir = os.path.join(protein_dir, "ligands")

        os.makedirs(ligands_dir, exist_ok=True)

        pdb_ids = group[PDB_COLUMN].dropna().astype(str).unique()

        receptor_saved = False

        for pdb_id in pdb_ids:
            pdb_id = pdb_id.strip().upper()

            pdb_file = download_pdb(pdb_id, protein_dir)
            if not pdb_file:
                continue

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, pdb_file)

            # ✅ Sauvegarde UN SEUL receptor
            if not receptor_saved:
                receptor_path = os.path.join(protein_dir, "receptor.pdb")
                save_receptor(structure, receptor_path)
                receptor_saved = True

            # ✅ Sauvegarde ligand spécifique
            ligand_path = os.path.join(ligands_dir, f"ligand_{pdb_id}.pdb")
            save_ligand(structure, ligand_path)

            # nettoyage ligand vide
            if os.path.getsize(ligand_path) < 300:
                os.remove(ligand_path)

    print("\n✅ Terminé !")
if __name__ == "__main__":
    main()
