#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/nettoyage_eau_%j.out
#SBATCH -e logs/nettoyage_eau_%j.err



cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome || exit 1

python3 - <<'PY'
import os

root = r"/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationEcoli_sansKiKd"

for target in sorted(os.listdir(root)):
    d = os.path.join(root, target)
    mol2_path = os.path.join(d, "crystal_ligand_clean.mol2")
    if not os.path.isfile(mol2_path):
        continue

    out_path = os.path.join(d, "crystal_ligand_clean.mol2")

    with open(mol2_path) as f:
        lines = f.readlines()

    atoms, bonds, headers = [], [], []
    in_atom = in_bond = False
    hoh_atom_ids = set()

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom, in_bond = True, False
            continue
        elif line.startswith("@<TRIPOS>BOND"):
            in_atom, in_bond = False, True
            continue
        elif line.startswith("@<TRIPOS>"):
            in_atom = in_bond = False
            headers.append(line)
            continue

        if in_atom:
            parts = line.split()
            if len(parts) >= 8:
                atom_id = int(parts[0])
                resname = parts[7]
                if "HOH" in resname or "WAT" in resname:
                    hoh_atom_ids.add(atom_id)
                    continue
                atoms.append(parts)
        elif in_bond:
            bonds.append(line.split())
        else:
            headers.append(line)

    # --- number ---
    idmap, new_atoms = {}, []
    for new_id, parts in enumerate(atoms, start=1):
        old_id = int(parts[0])
        idmap[old_id] = new_id
        parts[0] = f"{new_id:7d}"
        new_atoms.append(" ".join(parts) + "\n")

    new_bonds = []
    for bond_id, parts in enumerate(bonds, start=1):
        if len(parts) < 4:
            continue
        a1, a2 = int(parts[1]), int(parts[2])
        if a1 in idmap and a2 in idmap:
            parts[0] = f"{bond_id:7d}"
            parts[1] = f"{idmap[a1]:5d}"
            parts[2] = f"{idmap[a2]:5d}"
            new_bonds.append(" ".join(parts) + "\n")

    n_atoms, n_bonds = len(new_atoms), len(new_bonds)

    with open(out_path, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{os.path.join(d, 'ligand.pdb')}\n")
        f.write(f"{n_atoms:6d} {n_bonds:6d} 0 0 0\n")
        f.write("SMALL\nGASTEIGER\n\n")
        f.write("@<TRIPOS>ATOM\n")
        f.writelines(new_atoms)
        f.write("@<TRIPOS>BOND\n")
        f.writelines(new_bonds)

    print(f"{target}: {len(hoh_atom_ids)} HOH/WAT deleted, {n_atoms} atoms and {n_bonds} bonds kept.")


PY
