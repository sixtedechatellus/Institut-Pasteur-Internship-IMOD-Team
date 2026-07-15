from rdkit import Chem
import os, sys

ROOT = "data/ValidationEcoli_sansKiKd"  # chemin racine de ton dataset

def clean_smiles(smi):
    """
    Supprime tout après la première molécule principale et retire les fragments métalliques,
    renvoie un SMILES RDKit canonique si possible.
    """
    parts = smi.split("\t")[0].split()[0]   # garde juste la première colonne
    # garde uniquement la première molécule avant les métaux ou séparateurs '.'
    parts = parts.split(".")[0]
    try:
        mol = Chem.MolFromSmiles(parts)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

for sub in sorted(os.listdir(ROOT)):
    subdir = os.path.join(ROOT, sub)
    f_in = os.path.join(subdir, "actives_final.ism")
    f_out = os.path.join(subdir, "actives_final.ism")
    if not os.path.exists(f_in): 
        continue

    smiles = []
    with open(f_in) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            smi = clean_smiles(line)
            if smi:
                smiles.append(smi)
            else:
                print(f"⚠️  Ignoré (invalide) dans {sub}: {line}")

    if smiles:
        with open(f_out, "w") as fout:
            for s in smiles:
                fout.write(s + "\n")
        print(f"✅  {len(smiles)} SMILES valides écrits → {f_out}")
    else:
        print(f"❌  Aucun SMILES valide dans {sub}")
