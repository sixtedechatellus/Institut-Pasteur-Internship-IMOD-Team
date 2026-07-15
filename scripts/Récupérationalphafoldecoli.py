#!/usr/bin/env python3
import os
import time
import pandas as pd
import requests

INPUT_CSV = "Résultats/Résultats finaux propres/Escherichia_Coli_K12_protéome_complet.csv"
OUT_DIR = "Résultats/Résultats finaux propres/Alphafold"
os.makedirs(OUT_DIR, exist_ok=True)

BAD_LIGANDS = {
    "water", "sodium", "chloride", "magnesium", "calcium",
    "glycerol", "peg", "ethanol", "methanol"
}


def get_alphafold_pdb_url(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return None
        data = r.json()
        if not data:
            return None
        return data[0].get("pdbUrl")
    except Exception as e:
        print(f"AlphaFold API error for {uniprot_id}: {e}")
        return None


def download_file(url, out_path):
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and len(r.content) > 1000:
            with open(out_path, "wb") as f:
                f.write(r.content)
            return True
    except Exception as e:
        print(f"Download error: {e}")
    return False


def get_chembl_ligands(uniprot_id):
    ligands = []
    try:
        url = (
            f"https://www.ebi.ac.uk/chembl/api/data/target.json"
            f"?target_components__accession={uniprot_id}"
        )
        r = requests.get(url, timeout=20)
        data = r.json()
        if not data.get("targets"):
            return ligands

        target_id = data["targets"][0]["target_chembl_id"]
        act_url = (
            f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
            f"?target_chembl_id={target_id}&limit=1000"
        )
        r2 = requests.get(act_url, timeout=30)
        acts = r2.json().get("activities", [])

        for a in acts:
            smi = a.get("canonical_smiles")
            name = a.get("molecule_chembl_id")
            if smi and len(smi) >= 5:
                ligands.append((name, smi))

    except Exception as e:
        print(f"ChEMBL error for {uniprot_id}: {e}")

    return ligands


df = pd.read_csv(INPUT_CSV)
n_done = 0
n_skipped = 0
n_failed = 0

for _, row in df.iterrows():
    uid = str(row["uniprot_id"]).strip()

    if row["has_pdb"]:
        continue

    ligands = get_chembl_ligands(uid)
    valid = [(name, smi) for name, smi in ligands if smi]

    if not valid:
        n_skipped += 1
        continue

    af_url = get_alphafold_pdb_url(uid)
    if af_url is None:
        n_failed += 1
        continue

    protein_dir = os.path.join(OUT_DIR, uid)
    os.makedirs(protein_dir, exist_ok=True)

    pdb_path = os.path.join(protein_dir, "receptor.pdb")
    if not download_file(af_url, pdb_path):
        n_failed += 1
        continue

    ism_path = os.path.join(protein_dir, "ligands.ism")
    with open(ism_path, "w") as f:
        for name, smi in valid:
            f.write(f"{smi} {name}\n")

    n_done += 1
    time.sleep(0.3)

print(f"Done: {n_done} | Skipped: {n_skipped} | Failed: {n_failed}")


#%% Récupération des .pdb des ligands pour chaque dossier

import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles

ROOT = "Résultats/Résultats finaux propres/Alphafold"


def smiles_to_pdb(smiles, chembl_id, out_path):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                return False

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        except Exception:
            pass

        mol = Chem.RemoveHs(mol)
        rdmolfiles.MolToPDBFile(mol, out_path)
        return True

    except Exception as e:
        print(f"RDKit error for {chembl_id}: {e}")
        return False


def read_ism(ism_path):
    entries = []
    with open(ism_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], parts[1]))
    return entries


n_dirs = 0
n_ok = 0
n_failed = 0
n_skipped = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    ism_path = os.path.join(subdir_path, "ligands.ism")
    if not os.path.exists(ism_path):
        continue

    n_dirs += 1
    entries = read_ism(ism_path)

    if not entries:
        n_skipped += 1
        continue

    for smiles, chembl_id in entries:
        out_path = os.path.join(subdir_path, f"ligand_{chembl_id}.pdb")

        if os.path.exists(out_path) and os.path.getsize(out_path) > 100:
            n_skipped += 1
            continue

        if smiles_to_pdb(smiles, chembl_id, out_path):
            n_ok += 1
        else:
            n_failed += 1

print(f"Dirs: {n_dirs} | Generated: {n_ok} | Skipped: {n_skipped} | Failed: {n_failed}")


#%% Enlever les doublons du ligand.ism

import os

ROOT = "Résultats/Résultats finaux propres/Alphafold"

n_dirs = 0
n_removed = 0
n_kept = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    ism_path = os.path.join(subdir_path, "ligands.ism")
    if not os.path.exists(ism_path):
        continue

    n_dirs += 1
    seen_smiles = set()
    seen_chembl = set()
    cleaned_lines = []
    removed = 0

    with open(ism_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                cleaned_lines.append(line + "\n")
                continue

            smiles = parts[0]
            chembl_id = parts[1]

            if smiles in seen_smiles or chembl_id in seen_chembl:
                removed += 1
                continue

            seen_smiles.add(smiles)
            seen_chembl.add(chembl_id)
            cleaned_lines.append(line + "\n")

    with open(ism_path, "w") as f:
        f.writelines(cleaned_lines)

    n_kept += len(cleaned_lines)
    n_removed += removed

print(f"Dirs: {n_dirs} | Kept: {n_kept} | Removed: {n_removed}")


#%% Fpocket

import os
import subprocess
import shutil
import tempfile

ROOT = "DUD-E_fpocket"
count_success = 0
count_fail = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    receptor_path = os.path.join(subdir_path, "receptor.pdb")
    if not os.path.exists(receptor_path):
        continue

    work_dir = os.path.join(subdir_path, "fpocket_tmp")
    os.makedirs(work_dir, exist_ok=True)
    temp_receptor = os.path.join(work_dir, "target.pdb")
    shutil.copy(receptor_path, temp_receptor)

    try:
        cmd = ["fpocket", "-f", "target.pdb", "-m", "40", "-M", "150"]
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=120)

        pockets_dir = os.path.join(work_dir, "target_out", "pockets")
        pocket_0_file = os.path.join(pockets_dir, "pocket0_atm.pdb")

        if os.path.exists(pocket_0_file):
            final_pocket_path = os.path.join(subdir_path, "pocket_detected.pdb")
            shutil.copy(pocket_0_file, final_pocket_path)
            count_success += 1
        else:
            print(f"No pocket found for {subdir}")
            count_fail += 1

    except subprocess.TimeoutExpired:
        print(f"Timeout: {subdir}")
        count_fail += 1
    except FileNotFoundError:
        print("fpocket not found. Run: conda install -c conda-forge fpocket")
        count_fail += 1
        break
    except Exception as e:
        print(f"Error {subdir}: {e}")
        count_fail += 1

print(f"Done: {count_success} | Failed: {count_fail}")


#%% Nettoyer out de fpocket

import os
import subprocess
import shutil
import glob
import re


def natural_sort_key(s):
    match = re.search(r'pocket(\d+)', s)
    if match:
        return int(match.group(1))
    return 9999


ROOT = "DUD-E_fpocket"
NB_POCKETS_TO_KEEP = 5
count_success = 0
count_fail = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    receptor_path = os.path.join(subdir_path, "receptor.pdb")
    if not os.path.exists(receptor_path):
        continue

    work_dir = os.path.join(subdir_path, "fpocket_tmp")
    os.makedirs(work_dir, exist_ok=True)
    temp_receptor = os.path.join(work_dir, "target.pdb")
    shutil.copy(receptor_path, temp_receptor)

    try:
        cmd = ["fpocket", "-f", "target.pdb"]
        subprocess.run(cmd, cwd=work_dir, capture_output=True)

        out_folders = glob.glob(os.path.join(work_dir, "*_out"))
        if not out_folders:
            count_fail += 1
            shutil.rmtree(work_dir)
            continue

        pockets_src_dir = os.path.join(out_folders[0], "pockets")
        if not os.path.exists(pockets_src_dir):
            count_fail += 1
            shutil.rmtree(work_dir)
            continue

        all_files = [f for f in os.listdir(pockets_src_dir) if f.startswith("pocket") and f.endswith("_atm.pdb")]
        if not all_files:
            count_fail += 1
            shutil.rmtree(work_dir)
            continue

        all_files.sort(key=natural_sort_key)
        files_to_copy = all_files[:NB_POCKETS_TO_KEEP]

        for i, src_filename in enumerate(files_to_copy):
            src_path = os.path.join(pockets_src_dir, src_filename)
            dest_path = os.path.join(subdir_path, f"pocket_{i}.pdb")
            shutil.copy(src_path, dest_path)

        count_success += 1
        shutil.rmtree(work_dir)

    except Exception as e:
        print(f"Error {subdir}: {e}")
        count_fail += 1
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)

print(f"Done: {count_success} | Failed: {count_fail}")


#%% Création du pocket.lmdb

import os
import lmdb
import pickle
import numpy as np
from biopandas.pdb import PandasPdb

ROOT = "DUD-E_fpocket"
POCKET_INDEX = 0
count_success = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    pocket_file = os.path.join(subdir_path, f"pocket_{POCKET_INDEX}.pdb")
    if not os.path.exists(pocket_file):
        continue

    try:
        pp = PandasPdb().read_pdb(pocket_file)
        df = pp.df['ATOM']

        if df.empty:
            continue

        coords = df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
        atom_types = [str(a).strip() for a in df['atom_name'].tolist()]

        pocket_data = {
            'pocket': subdir,
            'pocket_index': POCKET_INDEX,
            "pocket_atoms": atom_types,
            "pocket_coordinates": coords
        }

        lmdb_path = os.path.join(subdir_path, "pocket.lmdb")
        env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False,
                        readahead=False, meminit=False, map_size=1099511627776)

        with env.begin(write=True) as txn:
            txn.put(str(0).encode('ascii'), pickle.dumps(pocket_data))

        env.close()
        count_success += 1

    except Exception as e:
        print(f"Error {subdir}: {e}")
        import traceback
        traceback.print_exc()

print(f"Done: {count_success} pocket.lmdb files created")


#%% Compter le nombre d'atomes dans lmdb

import lmdb
import pickle
import os

lmdb_path = "Résultats/Résultats finaux propres/Alphafold/P26647/pocket.lmdb"

env = lmdb.open(lmdb_path, subdir=False, readonly=True, lock=False)
with env.begin() as txn:
    data = pickle.loads(txn.get(b'0'))
    coords = data['pocket_coordinates']
    atoms = data['pocket_atoms']
    print(f"File size: {os.path.getsize(lmdb_path) / 1024:.1f} KB")
    print(f"Atoms in pocket: {len(coords)}")
    print(f"Atom types (sample): {atoms[:5]}")
env.close()


#%% Nettoyage fichiers

import os
import subprocess
import shutil
import glob

ROOT = "Résultats/Résultats finaux propres/Alphafold"

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    receptor_path = os.path.join(subdir_path, "receptor.pdb")
    if not os.path.exists(receptor_path):
        continue

    for f in glob.glob(os.path.join(subdir_path, "pocket_*.pdb")):
        os.remove(f)

    work_dir = os.path.join(subdir_path, "fpocket_tmp")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)

    os.makedirs(work_dir)
    temp_rec = os.path.join(work_dir, "target.pdb")
    shutil.copy(receptor_path, temp_rec)

    try:
        cmd = ["fpocket", "-f", "target.pdb"]
        subprocess.run(cmd, cwd=work_dir, capture_output=True, timeout=120)

        out_folders = glob.glob(os.path.join(work_dir, "*_out"))
        if out_folders:
            pockets_dir = os.path.join(out_folders[0], "pockets")
            if os.path.exists(pockets_dir):
                files = sorted([f for f in os.listdir(pockets_dir) if f.endswith("_atm.pdb")])
                for i, f in enumerate(files[:10]):
                    src = os.path.join(pockets_dir, f)
                    dst = os.path.join(subdir_path, f"pocket_{i}.pdb")
                    shutil.copy(src, dst)
                print(f"{subdir}: {len(files[:10])} pockets extracted")

        shutil.rmtree(work_dir)

    except Exception as e:
        print(f"Error {subdir}: {e}")


#%% Extension de la poche

import os
import lmdb
import pickle
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist

ROOT = "DUD-E_fpocket"
EXPANSION_RADIUS = 6.0
count_success = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    pocket_src_file = os.path.join(subdir_path, "pocket_0.pdb")
    receptor_file = os.path.join(subdir_path, "receptor.pdb")

    if not os.path.exists(pocket_src_file) or not os.path.exists(receptor_file):
        continue

    try:
        pp_pocket = PandasPdb().read_pdb(pocket_src_file)
        if pp_pocket.df['ATOM'].empty:
            continue
        coords_pocket = pp_pocket.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].to_numpy()

        pp_rec = PandasPdb().read_pdb(receptor_file)
        df_rec = pp_rec.df['ATOM']
        if df_rec.empty:
            continue

        coords_rec = df_rec[['x_coord', 'y_coord', 'z_coord']].to_numpy()
        distances = cdist(coords_rec, coords_pocket, metric='euclidean')
        min_dists = np.min(distances, axis=1)
        mask = min_dists <= EXPANSION_RADIUS
        expanded_df = df_rec[mask].copy()

        final_coords = expanded_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
        final_atoms = expanded_df['atom_name'].astype(str).str.strip().tolist()

        pocket_data = {
            'pocket': subdir,
            'pocket_index': 0,
            "pocket_atoms": final_atoms,
            "pocket_coordinates": final_coords
        }

        lmdb_path = os.path.join(subdir_path, "pocket.lmdb")
        if os.path.exists(lmdb_path):
            os.remove(lmdb_path)

        env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False,
                        readahead=False, meminit=False, map_size=1099511627776)

        with env.begin(write=True) as txn:
            txn.put(str(0).encode('ascii'), pickle.dumps(pocket_data))

        env.close()
        count_success += 1

    except Exception as e:
        print(f"Error {subdir}: {e}")
        import traceback
        traceback.print_exc()

print(f"Done: {count_success} expanded pockets created")


#%% Création du mols.lmdb

import os
import lmdb
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = "data/Alphafold"


def get_mol_data(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None

        mol = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        if res == -1:
            res = AllChem.EmbedMolecule(mol, randomSeed=42)
            if res == -1:
                return None

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except:
            pass

        mol = Chem.RemoveHs(mol)

        if mol.GetNumConformers() == 0:
            return None

        coords = mol.GetConformer().GetPositions()
        atom_types = [a.GetSymbol() for a in mol.GetAtoms()]

        return {
            'atoms': atom_types,
            'coordinates': [coords],
            'smi': Chem.MolToSmiles(mol),
        }
    except Exception:
        return None


count_success = 0
total_actives = 0
total_decoys = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    actives_file = os.path.join(subdir_path, "actives_final.ism")
    decoys_file = os.path.join(subdir_path, "decoys_final.ism")

    if not os.path.exists(actives_file):
        continue

    all_mol_entries = []

    n_act = 0
    with open(actives_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                smi = parts[0]
                name = parts[1] if len(parts) > 1 else "Active"
                data = get_mol_data(smi)
                if data:
                    data['label'] = 1
                    data['name'] = name
                    all_mol_entries.append(data)
                    n_act += 1

    n_dec = 0
    if os.path.exists(decoys_file):
        with open(decoys_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    smi = parts[1]
                    name = f"Decoy_{n_dec}"
                elif len(parts) == 1:
                    smi = parts[0]
                    name = f"Decoy_{n_dec}"
                else:
                    continue
                data = get_mol_data(smi)
                if data:
                    data['label'] = 0
                    data['name'] = name
                    all_mol_entries.append(data)
                    n_dec += 1

    if not all_mol_entries:
        continue

    lmdb_path = os.path.join(subdir_path, "mols.lmdb")
    if os.path.exists(lmdb_path):
        os.remove(lmdb_path)

    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False,
                    readahead=False, meminit=False, map_size=1099511627776)

    with env.begin(write=True) as txn:
        for idx, entry in enumerate(all_mol_entries):
            txn.put(str(idx).encode('ascii'), pickle.dumps(entry))

    env.close()

    total_actives += n_act
    total_decoys += n_dec
    count_success += 1

print(f"Done: {count_success} dirs | Actives: {total_actives} | Decoys: {total_decoys}")


#%% Nettoyage fichier smiles

import os

ROOT = "Résultats/Résultats finaux propres/Alphafold"
count_success = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    input_file = os.path.join(subdir_path, "ligands.ism")
    output_file = os.path.join(subdir_path, "actives_final.ism")

    if not os.path.exists(input_file):
        continue

    cleaned_lines = []

    try:
        with open(input_file, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                smiles = parts[0]
                if len(smiles) > 3:
                    cleaned_lines.append(smiles)

        if cleaned_lines:
            with open(output_file, 'w') as f_out:
                for smi in cleaned_lines:
                    f_out.write(smi + "\n")
            count_success += 1

    except Exception as e:
        print(f"Error {subdir}: {e}")

print(f"Done: {count_success} actives_final.ism files created")


#%% Poches P2rank

#!/usr/bin/env python3
import os
import subprocess
import shutil
import glob
import gzip
import pickle
import numpy as np
import lmdb
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist

ROOT = "data/Alphafold_P2rank"
EXPANSION_RADIUS = 6.0
PRANK_CMD = "/home/sixte/Documents/Helix/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/p2rank/prank"


def extract_first_model(in_file, out_file):
    lines = open(in_file).read().splitlines()
    with open(out_file, "w") as out:
        write = False
        for line in lines:
            if line.startswith("MODEL") and "1" in line.split():
                write = True
                continue
            if line.startswith("ENDMDL") and write:
                break
            if write:
                out.write(line + "\n")


count_success = 0

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    receptor_path = os.path.join(subdir_path, "receptor.pdb")
    if not os.path.exists(receptor_path):
        continue

    work_dir = os.path.join(subdir_path, "p2rank_tmp")
    os.makedirs(work_dir, exist_ok=True)
    temp_pdb = os.path.join(work_dir, "target.pdb")
    shutil.copy(receptor_path, temp_pdb)

    try:
        cmd = [
            PRANK_CMD, "predict",
            "-f", "target.pdb",
            "-o", "output", ".",
            "-visualizations", "1"
        ]
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, timeout=240)

        if result.returncode != 0:
            print(f"P2Rank error for {subdir}:\n{result.stderr}")
            shutil.rmtree(work_dir)
            continue

        vis_dir = os.path.join(work_dir, "output", "visualizations", "data")
        gz_file = glob.glob(os.path.join(vis_dir, "*_points.pdb.gz"))
        if not gz_file:
            shutil.rmtree(work_dir)
            continue

        gz_file = gz_file[0]
        decompressed = os.path.join(work_dir, "target_points.pdb")
        with gzip.open(gz_file, "rb") as f_in, open(decompressed, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        pocket_pdb = os.path.join(subdir_path, "pocket_0.pdb")
        shutil.copy(decompressed, pocket_pdb)

        pp_pocket = PandasPdb().read_pdb(pocket_pdb)
        df_pocket = pp_pocket.df.get("ATOM")
        if df_pocket is None or df_pocket.empty:
            df_pocket = pp_pocket.df.get("HETATM")
        if df_pocket is None or df_pocket.empty:
            shutil.rmtree(work_dir)
            continue

        coords_pocket = df_pocket[["x_coord", "y_coord", "z_coord"]].to_numpy()

        pp_rec = PandasPdb().read_pdb(receptor_path)
        df_rec = pp_rec.df.get("ATOM")
        if df_rec is None or df_rec.empty:
            df_rec = pp_rec.df.get("HETATM")
        if df_rec is None or df_rec.empty:
            rx, ry, rz, rnames = [], [], [], []
            with open(receptor_path) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                        try:
                            rx.append(float(line[30:38]))
                            ry.append(float(line[38:46]))
                            rz.append(float(line[46:54]))
                            rnames.append(line[12:16].strip())
                        except ValueError:
                            continue
            df_rec = pd.DataFrame({
                "x_coord": rx, "y_coord": ry, "z_coord": rz, "atom_name": rnames
            })

        if df_rec.empty:
            shutil.rmtree(work_dir)
            continue

        coords_rec = df_rec[["x_coord", "y_coord", "z_coord"]].to_numpy()
        distances = cdist(coords_rec, coords_pocket)
        min_dists = np.min(distances, axis=1)
        mask = min_dists <= EXPANSION_RADIUS
        expanded_df = df_rec[mask].copy()

        final_coords = expanded_df[["x_coord", "y_coord", "z_coord"]].to_numpy()
        final_atoms = expanded_df["atom_name"].astype(str).str.strip().tolist()

        pocket_data = {
            "pocket": subdir,
            "pocket_index": 0,
            "pocket_atoms": final_atoms,
            "pocket_coordinates": final_coords,
        }

        lmdb_path = os.path.join(subdir_path, "pocket.lmdb")
        if os.path.exists(lmdb_path):
            os.remove(lmdb_path)

        env = lmdb.open(lmdb_path, subdir=False, lock=False, map_size=1099511627776)
        with env.begin(write=True) as txn:
            txn.put(b"0", pickle.dumps(pocket_data))
        env.close()

        count_success += 1
        shutil.rmtree(work_dir)

    except Exception as e:
        print(f"Error {subdir}: {e}")
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)

print(f"Done: {count_success} pockets generated")


#%% Fonction nécessaire pour script précédent

import pandas as pd


def read_points_pdb(file_path):
    """Lit un fichier PDB de points P2Rank, ignore les molécules d'eau."""
    atoms = []
    with open(file_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
                resname = line[17:20].strip().upper()
                if resname in {"HOH", "WAT", "H2O"}:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append((x, y, z))
                except ValueError:
                    continue

    df = pd.DataFrame(atoms, columns=["x_coord", "y_coord", "z_coord"])
    return df


#%% Vérification pocket.lmdb

import lmdb
import pickle
import numpy as np
import os

lmdb_path = "DUD-E/aa2ar/pocket.lmdb"

if not os.path.exists(lmdb_path):
    print(f"File not found: {lmdb_path}")
    exit(1)

env = lmdb.open(lmdb_path, subdir=False, readonly=True, lock=False, readahead=False, meminit=False)
with env.begin() as txn:
    stats = env.stat()
    print(f"Entries: {stats['entries']}")

    cursor = txn.cursor()
    for key, value in cursor:
        data = pickle.loads(value)
        print(f"Keys: {list(data.keys())}")

        if 'pocket_atoms' in data:
            atoms = data['pocket_atoms']
            print(f"pocket_atoms: {len(atoms)} atoms, sample: {atoms[:5]}")
        else:
            print("Missing key: pocket_atoms")

        if 'pocket_coordinates' in data:
            coords = data['pocket_coordinates']
            shape = coords.shape if hasattr(coords, 'shape') else 'N/A'
            print(f"pocket_coordinates: shape={shape}")
        else:
            print("Missing key: pocket_coordinates")

        break

env.close()


#%% Shuffle des decoys

import os
import random
import shutil

main_folder = r"data/ValidationEcoli_sansKiKd_shuffle"
decoys_filename = "decoys_final.ism"
create_backup = True
random_seed = 42


def read_decoys_file(filepath):
    pairs = []
    skipped = 0
    with open(filepath, 'r') as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split(' ')
            if len(parts) == 2:
                pairs.append((parts[0], parts[1]))
            else:
                skipped += 1
    return pairs, skipped


def write_decoys_file(filepath, pairs):
    with open(filepath, 'w') as f:
        for actif, decoy in pairs:
            f.write(f"{actif} {decoy}\n")


if not os.path.exists(main_folder):
    raise FileNotFoundError(f"Folder not found: {main_folder}")

subdirs = sorted([
    d for d in os.listdir(main_folder)
    if os.path.isdir(os.path.join(main_folder, d))
])

all_files_pairs = {}
all_decoys = []

for subdir in subdirs:
    decoys_file = os.path.join(main_folder, subdir, decoys_filename)

    if not os.path.exists(decoys_file):
        continue

    if create_backup:
        backup_path = decoys_file + ".backup"
        if not os.path.exists(backup_path):
            shutil.copy2(decoys_file, backup_path)

    pairs, skipped = read_decoys_file(decoys_file)
    if not pairs:
        continue

    all_files_pairs[decoys_file] = pairs
    all_decoys.extend(decoy for _, decoy in pairs)

if len(all_decoys) == 0:
    raise ValueError("No decoys found.")

if random_seed is not None:
    random.seed(random_seed)

random.shuffle(all_decoys)

decoy_idx = 0
for filepath, pairs in all_files_pairs.items():
    new_pairs = []
    for actif, _ in pairs:
        new_pairs.append((actif, all_decoys[decoy_idx]))
        decoy_idx += 1
    write_decoys_file(filepath, new_pairs)

print(f"Done: {len(all_files_pairs)} files | {decoy_idx} decoys redistributed")


#%% Filtre enlève les decoys qui ont 5 atomes plus lourds que la molécule active

import os
from rdkit import Chem

main_folder = r"data/ValidationEcoli_sansKiKd_shuffle"
decoys_filename = "decoys_final.ism"
max_heavy_diff = 5


def count_heavy_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol.GetNumHeavyAtoms()


def filter_file(filepath, max_diff):
    pairs_ok = []
    nb_supprime = 0
    nb_invalide = 0

    with open(filepath, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    for line in lines:
        parts = line.split(' ')
        if len(parts) != 2:
            nb_invalide += 1
            continue

        actif, decoy = parts
        n_actif = count_heavy_atoms(actif)
        n_decoy = count_heavy_atoms(decoy)

        if n_actif is None or n_decoy is None:
            nb_invalide += 1
            continue

        if abs(n_actif - n_decoy) <= max_diff:
            pairs_ok.append((actif, decoy))
        else:
            nb_supprime += 1

    with open(filepath, 'w') as f:
        for actif, decoy in pairs_ok:
            f.write(f"{actif} {decoy}\n")

    return len(lines) - nb_invalide, len(pairs_ok), nb_supprime, nb_invalide


if not os.path.exists(main_folder):
    raise FileNotFoundError(f"Folder not found: {main_folder}")

subdirs = sorted([d for d in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, d))])

total_initial = 0
total_garde = 0
total_supprime = 0
total_invalide = 0
fichiers_traites = 0

for subdir in subdirs:
    filepath = os.path.join(main_folder, subdir, decoys_filename)
    if not os.path.exists(filepath):
        continue

    nb_init, nb_ok, nb_sup, nb_inv = filter_file(filepath, max_heavy_diff)
    total_initial += nb_init
    total_garde += nb_ok
    total_supprime += nb_sup
    total_invalide += nb_inv
    fichiers_traites += 1

print(f"Files: {fichiers_traites}")
print(f"Initial: {total_initial} | Kept: {total_garde} | Removed: {total_supprime} | Invalid: {total_invalide}")
if total_initial > 0:
    print(f"Retention rate: {total_garde / total_initial * 100:.1f}%")


#%% Nettoyer les dossiers sans decoys

import os
import shutil

main_folder = r"data/ValidationEcoli_sansKiKd_shuffle"
target_filename = "decoys_final.ism"
dry_run = False

if not os.path.exists(main_folder):
    raise FileNotFoundError(f"Folder not found: {main_folder}")

subdirs = sorted([
    d for d in os.listdir(main_folder)
    if os.path.isdir(os.path.join(main_folder, d))
])

count_to_delete = 0
count_kept = 0

for subdir in subdirs:
    subdir_path = os.path.join(main_folder, subdir)
    target_file = os.path.join(subdir_path, target_filename)

    if os.path.exists(target_file):
        count_kept += 1
    else:
        count_to_delete += 1
        if not dry_run:
            try:
                shutil.rmtree(subdir_path)
            except Exception as e:
                print(f"Error deleting {subdir}: {e}")

print(f"Kept: {count_kept} | Deleted: {count_to_delete}")
if dry_run and count_to_delete > 0:
    print(f"Set dry_run=False to actually delete {count_to_delete} folders.")