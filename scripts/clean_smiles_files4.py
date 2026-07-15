from rdkit import Chem
import os

# ===============================
# CONFIG
# ===============================
ROOT = "data/ValidationEcoli_sansKiKd"
MIN_HEAVY_ATOMS = 6

ALLOWED_ATOMIC_SYMBOLS = {
    'C', 'H', 'N', 'O', 'S',
    'F', 'Cl', 'Br', 'I'
}

# ===============================
# FONCTION DE NETTOYAGE
# ===============================
def clean_smiles_for_deepcoy(raw_line):
    """
    Nettoie le SMILES pour DeepCOY :
    1. Split manuel par '.' avant RDKit (gère les complexes Fe-S illisibles).
    2. Lit chaque fragment indépendamment.
    3. Filtre taille, atomes interdits, solvants.
    4. Double canonicalisation (corrige [PH] → P standard).
    5. Retourne le plus gros fragment valide.
    """
    try:
        raw_smi = raw_line.split("\t")[0].split()[0].strip()
        if not raw_smi:
            return None

        # Split manuel par '.' AVANT RDKit
        raw_parts = raw_smi.split('.')
        candidates = []

        for part in raw_parts:
            if not part:
                continue

            # Lecture individuelle de chaque fragment
            frag = Chem.MolFromSmiles(part)
            if frag is None:
                continue  # Fragment illisible (ex: cluster Fe-S) → ignoré

            # 1. Filtre Taille
            if frag.GetNumHeavyAtoms() < MIN_HEAVY_ATOMS:
                continue

            # 2. Filtre Atomes (LISTE BLANCHE)
            is_valid = True
            for atom in frag.GetAtoms():
                if atom.GetSymbol() not in ALLOWED_ATOMIC_SYMBOLS:
                    is_valid = False
                    break
            if not is_valid:
                continue

            # 3. Double canonicalisation (corrige [PH] → P standard)
            temp_smi = Chem.MolToSmiles(frag, isomericSmiles=True)
            temp_mol = Chem.MolFromSmiles(temp_smi)
            if temp_mol is None:
                continue
            smi_frag = Chem.MolToSmiles(temp_mol, isomericSmiles=True)
            # Force remplacement [PH] → P si encore présent
            if '[PH]' in smi_frag:
                smi_frag = smi_frag.replace('[PH]', 'P')
                check = Chem.MolFromSmiles(smi_frag)
                if check is None:
                    continue
            # 4. Filtre Solvants
            if smi_frag in ["OCC(O)CO", "CS(=O)C", "CCO", "CO", "CC(C)=O"]:
                continue

            candidates.append((frag.GetNumHeavyAtoms(), smi_frag))

        if not candidates:
            return None

        # 5. Garder le plus gros fragment valide
        candidates.sort(key=lambda x: x[0], reverse=True)
        final_smi = candidates[0][1]

        # 6. Vérification finale
        if Chem.MolFromSmiles(final_smi) is None:
            return None

        return final_smi

    except Exception:
        return None



# ===============================
# MAIN LOOP
# ===============================
total_processed = 0
total_cleaned = 0
total_rejected = 0


for sub in sorted(os.listdir(ROOT)):
    subdir = os.path.join(ROOT, sub)
    if not os.path.isdir(subdir):
        continue

    input_file = os.path.join(subdir, "actives_final.ism")
    output_file = input_file

    if not os.path.exists(input_file):
        continue

    unique_smiles = set()
    local_rejected = 0

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            cleaned = clean_smiles_for_deepcoy(line)

            if cleaned:
                unique_smiles.add(cleaned)
            else:
                local_rejected += 1

    if unique_smiles:

        with open(output_file, 'w') as f:
            for smi in sorted(list(unique_smiles)):
                f.write(smi + "\n")

        kept = len(unique_smiles)
        total_cleaned += kept
        total_rejected += local_rejected
    else:
        if os.path.exists(output_file):
            os.remove(output_file)
        total_rejected += local_rejected

    total_processed += 1

print(f"Dossiers traités     : {total_processed}")
print(f"Ligands propres      : {total_cleaned}")
print(f"Ligands rejetés      : {total_rejected}")

