import os
import pandas as pd
import requests
import re
from tqdm import tqdm

# ===============================
# CONFIG
# ===============================
CSV_PATH = "/pasteur/helix/projects/Imod-moulin/Sixte/Code/Récupération_base_de_données_Chembl/Résultats/Résultats_finaux_propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Biolip2.csv"
OUTPUT_DIR = "data/ValidationEcoli_sansKiKd"

PDB_COLUMN = "PDB_ID"
PROTEIN_COLUMN = "UniProt_ID"

# ===============================
# DOWNLOAD & PROCESS (Optimisé sans BioPython lourd)
# ===============================

def download_and_process_pdb(pdb_id, protein_dir, ligands_dir, receptor_saved):
    """
    Télécharge le PDB et extrait directement les lignes ATOM et HETATM via Regex.
    Beaucoup plus rapide que Bio.PDB pour ce cas d'usage.
    """
    pdb_id = pdb_id.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            return False, receptor_saved
        
        content = r.text
        if "ATOM" not in content:
            return False, receptor_saved

        lines = content.splitlines()
        
        # 1. Sauvegarde du Récepteur (uniquement si pas encore fait pour cette protéine)
        if not receptor_saved:
            receptor_path = os.path.join(protein_dir, "receptor.pdb")
            with open(receptor_path, "w") as f:
                for line in lines:
                    if line.startswith("ATOM"):
                        f.write(line + "\n")
            receptor_saved = True

        # 2. Sauvegarde du Ligand spécifique à ce PDB
        ligand_lines = []
        for line in lines:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                # Ignorer eau et sels communs
                if resname not in ["HOH", "WAT", "Na", "Cl", "Mg", "Ca", "K"]:
                    ligand_lines.append(line)
        
        if ligand_lines:
            ligand_path = os.path.join(ligands_dir, f"ligand_{pdb_id}.pdb")
            with open(ligand_path, "w") as f:
                # On ajoute un header minimal pour que ce soit un vrai PDB
                f.write(f"HEADER    LIGAND FOR {pdb_id}\n")
                f.write("\n".join(ligand_lines) + "\n")
                f.write("END\n")
            
            # Nettoyage si trop petit (moins de 5 lignes)
            if len(ligand_lines) < 5:
                os.remove(ligand_path)
                
        return True, receptor_saved

    except Exception as e:
        # print(f"Erreur sur {pdb_id}: {e}")
        return False, receptor_saved

# ===============================
# MAIN
# ===============================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Correction du warning : on force la colonne 16 en string si nécessaire, ou on ignore le warning
    # Souvent la colonne 'Ligand_Name' ou similaire contient des chiffres et lettres
    try:
        df = pd.read_csv(CSV_PATH, dtype=str, low_memory=False)
    except Exception as e:
        print(f"Erreur lecture CSV: {e}")
        return

    # Vérification des colonnes
    if PDB_COLUMN not in df.columns or PROTEIN_COLUMN not in df.columns:
        print(f"Colonnes manquantes. Disponibles: {df.columns}")
        return

    # Nettoyage des IDs
    df[PDB_COLUMN] = df[PDB_COLUMN].astype(str).str.strip().str.upper()
    df = df[df[PDB_COLUMN] != "NAN"] # Retirer les NaN convertis en string

    grouped = df.groupby(PROTEIN_COLUMN)
    total_proteins = len(grouped)

    print(f"🔍 {total_proteins} protéines uniques à traiter")

    success_count = 0

    # Boucle principale
    for protein_id, group in tqdm(grouped, desc="Protéines", unit="prot"):
        
        protein_dir = os.path.join(OUTPUT_DIR, str(protein_id))
        ligands_dir = os.path.join(protein_dir, "ligands")
        os.makedirs(ligands_dir, exist_ok=True)

        pdb_ids = group[PDB_COLUMN].unique()
        receptor_saved = False

        for pdb_id in pdb_ids:
            if not pdb_id or len(pdb_id) != 4:
                continue
                
            ok, receptor_saved = download_and_process_pdb(
                pdb_id, protein_dir, ligands_dir, receptor_saved
            )
            if ok:
                success_count += 1

    print(f"\n✅ Terminé ! {success_count} fichiers PDB traités avec succès.")

if __name__ == "__main__":
    main()
