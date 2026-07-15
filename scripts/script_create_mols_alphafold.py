#!/usr/bin/env python3
"""
Crée les fichiers mols.lmdb individuels pour chaque protéine.
Optimisé pour exécution sur cluster (pas de barre de progression interactive lourde).
"""

import os
import sys
import lmdb
import pickle
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# Désactiver les warnings RDKit (UFFTYPER, etc.) pour des logs propres
RDLogger.DisableLog('rdApp.*')

# ============================================================
# ⚙️ CONFIG
# ============================================================
# Prend le dossier en argument ou utilise une valeur par défaut
if len(sys.argv) > 1:
    ROOT = sys.argv[1]
else:
    ROOT = "data/Alphafold"  # Change ceci si besoin

print(f"🚀 Démarrage création mols.lmdb dans : {ROOT}")

# ============================================================
# FONCTIONS RDKit
# ============================================================
def get_mol_data(smiles):
    """Génère les données 3D d'une molécule depuis son SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        mol = Chem.AddHs(mol)
        # ETKDGv3 est robuste
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

# ============================================================
# MAIN
# ============================================================
count_success = 0
total_actives = 0
total_decoys = 0
failed_dirs = 0

if not os.path.exists(ROOT):
    print(f"❌ Erreur : Le dossier {ROOT} n'existe pas.")
    sys.exit(1)

for subdir in sorted(os.listdir(ROOT)):
    subdir_path = os.path.join(ROOT, subdir)
    if not os.path.isdir(subdir_path):
        continue

    actives_file = os.path.join(subdir_path, "actives_final.ism") 
    decoys_file = os.path.join(subdir_path, "decoys_final.ism")
    
    if not os.path.exists(actives_file):
        continue

    print(f"[→] Traitement : {subdir}...", flush=True)

    all_mol_entries = []
    n_act = 0
    n_dec = 0

    # 1. Charger les ACTIFS (Label = 1)
    try:
        with open(actives_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 1:
                    smi = parts[0]
                    data = get_mol_data(smi)
                    if data:
                        data['label'] = 1
                        data['name'] = parts[1] if len(parts) > 1 else "Active"
                        all_mol_entries.append(data)
                        n_act += 1
    except Exception as e:
        print(f"   ⚠️ Erreur lecture actifs: {e}")

    # 2. Charger les DECOYS (Label = 0)
    # Gestion intelligente : si 2 colonnes, prend la 2ème. Si 1, prend la 1ère.
    if os.path.exists(decoys_file):
        try:
            with open(decoys_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    smi = None
                    if len(parts) >= 2:
                        smi = parts[1] # Cas DeepCoy : 2ème colonne est le décoy
                    elif len(parts) == 1:
                        smi = parts[0] # Cas standard
                    
                    if smi:
                        data = get_mol_data(smi)
                        if data:
                            data['label'] = 0
                            data['name'] = f"Decoy_{n_dec}"
                            all_mol_entries.append(data)
                            n_dec += 1
        except Exception as e:
            print(f"   ⚠️ Erreur lecture décoys: {e}")

    if not all_mol_entries:
        print(f"   ⚠️ Aucune molécule valide. Skip.", flush=True)
        failed_dirs += 1
        continue

    # 3. Écrire dans mols.lmdb
    lmdb_path = os.path.join(subdir_path, "mols.lmdb")
    if os.path.exists(lmdb_path):
        os.remove(lmdb_path)
    
    try:
        env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False, 
                        readahead=False, meminit=False, map_size=1099511627776)
        
        with env.begin(write=True) as txn:
            for idx, entry in enumerate(all_mol_entries):
                txn.put(str(idx).encode('ascii'), pickle.dumps(entry))
        env.close()
        
        print(f"   ✅ OK : {n_act} actifs, {n_dec} décoys.", flush=True)
        count_success += 1
        total_actives += n_act
        total_decoys += n_dec
        
    except Exception as e:
        print(f"   ❌ Erreur écriture LMDB: {e}", flush=True)
        failed_dirs += 1

print(f"\n{'='*50}")
print(f"TERMINÉ !")
print(f"Dossiers succès : {count_success}")
print(f"Dossiers échec/vides : {failed_dirs}")
print(f"Total Actifs traités : {total_actives}")
print(f"Total Décoys traités : {total_decoys}")
