import lmdb
import pickle
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

data_dir = "./data/lit_pcba/lit_pcba"
output_path = os.path.join(data_dir, "mols.lmdb")

# Read SMILES
smi_files = [f for f in os.listdir(data_dir) if f.endswith(".smi")]
all_mols = []

for smi_file in smi_files:
    with open(os.path.join(data_dir, smi_file)) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 1:
                smi = parts[0]
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)
                res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                if res != 0:
                    continue
                AllChem.MMFFOptimizeMolecule(mol)
                atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
                coords = mol.GetConformer().GetPositions().tolist()
                all_mols.append({
                    "smi": smi,
                    "atoms": atoms,
                    "coordinates": [coords],
                })

# write LMDB
env = lmdb.open(output_path, map_size=int(1e12))
with env.begin(write=True) as txn:
    for i, mol_data in enumerate(all_mols):
        txn.put(str(i).encode(), pickle.dumps(mol_data))
    txn.put(b"__len__", pickle.dumps(len(all_mols)))

print(f"✅ {len(all_mols)} molecules encoded in {output_path}")
