#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/nettoyage_eau_%j.out
#SBATCH -e logs/nettoyage_eau_%j.err

ROOT="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationEcoli_sansKiKd"

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome || exit 1

python - <<PY
import os, re

root = r"$ROOT"

for target in sorted(os.listdir(root)):
    d = os.path.join(root, target)
    mol2_path = os.path.join(d, "crystal_ligand.mol2")
    if not os.path.isfile(mol2_path):
        continue

    out_path = os.path.join(d, "crystal_ligand_clean.mol2")

    with open(mol2_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    hoh_atom_ids = set()
    in_atom = False
    in_bond = False

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            in_bond = False
            new_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>BOND"):
            in_atom = False
            in_bond = True
            new_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            in_bond = False
            new_lines.append(line)
            continue

        if in_atom:
            # champs : ID, nom, x, y, z, type, resid, resname, charge
            parts = line.split()
            if len(parts) >= 8:
                atom_id = int(parts[0])
                resname = parts[7]
                if "HOH" in resname or "WAT" in resname:
                    hoh_atom_ids.add(int(atom_id))
                    # on n'ajoute pas cette ligne
                    continue
            new_lines.append(line)
        elif in_bond:
            # format : bond_id  a1  a2  type
            parts = line.split()
            if len(parts) >= 4:
                a1, a2 = int(parts[1]), int(parts[2])
                if a1 in hoh_atom_ids or a2 in hoh_atom_ids:
                    # bond lié à un atome supprimé -> on le retire
                    continue
            new_lines.append(line)
        else:
            # toutes les autres lignes (MOLECULE, CHARGE, etc.)
            new_lines.append(line)

    if hoh_atom_ids:
        with open(out_path, "w") as f:
            f.writelines(new_lines)
        print(f"[OK] {target}: supprimé {len(hoh_atom_ids)} atomes HOH/WAT et leurs bonds -> crystal_ligand_clean.mol2")
    else:
        print(f"[INFO] {target}: aucun HOH/WAT trouvé, rien changé.")

PY

echo "[✔] Nettoyage terminé."
