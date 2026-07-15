import os
import pandas as pd
import requests
from tqdm import tqdm
from Bio.PDB import PDBParser, PDBIO, Select

# ===============================
# CONFIGURATION
# ===============================
CSV_PATH = "/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/Validationdatasets/BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv"  # <-- ton fichier CSV
OUTPUT_DIR = "data/ValidationEcoli"  # dossier parent pour sauvegarder
PDB_COLUMN = "PDB_ID"             # nom de la colonne dans ton CSV

# ===============================
# UTILS
# ===============================
class ProteinSelect(Select):
    """Sélectionne uniquement les chaînes protéiques (ATOM)."""
    def accept_residue(self, residue):
        return residue.id[0] == " "  # résidus standards

class LigandSelect(Select):
    """Sélectionne les résidus non protéiques (ligands potentiels)."""
    def accept_residue(self, residue):
        if residue.id[0].startswith("H"):
            if residue.resname not in ["HOH", "WAT"]:
                return True
        return False

def download_pdb(pdb_id, dest_folder):
    """Télécharge un fichier PDB depuis RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    dest_file = os.path.join(dest_folder, f"{pdb_id}.pdb")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and "ATOM" in r.text:
            with open(dest_file, "w") as f:
                f.write(r.text)
            return dest_file
        else:
            print(f"⚠️  Fichier {pdb_id} non trouvé (status {r.status_code})")
            return None
    except Exception as e:
        print(f"Erreur téléchargement {pdb_id}: {e}")
        return None

def split_protein_ligand(pdb_path, out_dir):
    """Extrait protéine et ligand du PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)

    base_id = os.path.splitext(os.path.basename(pdb_path))[0]
    io = PDBIO()
    io.set_structure(structure)

    # récepteur (ATOM)
    receptor_path = os.path.join(out_dir, f"receptor_{base_id}.pdb")
    io.save(receptor_path, ProteinSelect())

    # ligands (HETATM sauf eau)
    ligand_path = os.path.join(out_dir, f"ligand_{base_id}.pdb")
    io.save(ligand_path, LigandSelect())

    return receptor_path, ligand_path

# ===============================
# MAIN logic
# ===============================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)
    if PDB_COLUMN not in df.columns:
        raise KeyError(f"Colonne '{PDB_COLUMN}' introuvable dans {CSV_PATH}")
    
    pdb_ids = df[PDB_COLUMN].dropna().astype(str).unique()

    print(f"🔍 {len(pdb_ids)} PDB IDs trouvés dans ton CSV.")
    for pdb_id in tqdm(pdb_ids, desc="Traitement des complexes"):
        pdb_id = pdb_id.strip().upper()
        
        # Crée un sous-dossier spécifique pour chaque PDB
        pdb_dir = os.path.join(OUTPUT_DIR, pdb_id)
        os.makedirs(pdb_dir, exist_ok=True)
        
        pdb_file = download_pdb(pdb_id, pdb_dir)
        if not pdb_file:
            continue

        try:
            receptor, ligand = split_protein_ligand(pdb_file, pdb_dir)
            # Vérifie que le ligand a du contenu
            if os.path.getsize(ligand) < 300:
                os.remove(ligand)
                print(f"⚠️  Aucun ligand clair dans {pdb_id}, supprimé.")
        except Exception as e:
            print(f"Erreur lors du parsing {pdb_id}: {e}")

    print("\n✅ Terminé !")
    print(f"Les fichiers générés se trouvent dans des dossiers séparés sous : {OUTPUT_DIR}/")
    print("Chaque dossier contient :")
    print("- receptor_<PDB_ID>.pdb")
    print("- ligand_<PDB_ID>.pdb")

if __name__ == "__main__":
    main()
