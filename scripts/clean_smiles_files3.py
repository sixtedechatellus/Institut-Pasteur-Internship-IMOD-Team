from rdkit import Chem
import os

# ===============================
# CONFIG
# ===============================
ROOT = "data/ValidationEcoli_sansKiKd"  # <-- Vérifie ton chemin
MIN_HEAVY_ATOMS = 6      # Taille minimale (évite fragments trop petits)

# 🛡️ LISTE BLANCHE DES ATOMES AUTORISÉS POUR DEEPCOY
# Seulement les atomes organiques standards trouvés dans ZINC/ChEMBL.
# Si un ligand contient un autre atome (B, Si, As, Se, métaux...), il sera rejeté.
ALLOWED_ATOMIC_SYMBOLS = {
    'C', 'H', 'N', 'O', 'S', 'P',  # Organique de base
    'F', 'Cl', 'Br', 'I'           # Halogènes standards
    # 'B', 'Si'  # Décommente uniquement si tu es SÛR que ton modèle DeepCOY les gère
}

# ===============================
# FONCTION DE NETTOYAGE
# ===============================
def clean_smiles_for_deepcoy(raw_line):
    """
    Nettoie le SMILES pour DeepCOY :
    1. Extrait et sépare les sels.
    2. Garde le plus gros fragment organique.
    3. Vérifie la taille minimale.
    4. Vérifie que TOUS les atomes sont dans la liste blanche.
    5. Rejette les solvants connus (Glycérol, etc.).
    6. Canonicalise.
    """
    try:
        raw_smi = raw_line.split("\t")[0].split()[0].strip()
        if not raw_smi:
            return None

        # On split d'abord manuellement par '.'
        # Cela permet de gérer les fragments illisibles (ex: complexes Fe-S)
        # sans rejeter toute la molécule
        raw_parts = raw_smi.split('.')
        
        candidates = []
        for part in raw_parts:
            if not part:
                continue
            # On essaie de lire chaque fragment indépendamment
            frag = Chem.MolFromSmiles(part)
            if frag is None:
                continue  # Fragment illisible (ex: Fe-S cluster) → on ignore juste ce fragment
        
        for frag in frags:
            num_heavy = frag.GetNumHeavyAtoms()
            
            # 1. Filtre Taille
            if num_heavy < MIN_HEAVY_ATOMS:
                continue
            
            # 2. Filtre Atomes (LISTE BLANCHE)
            is_valid_composition = True
            for atom in frag.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol not in ALLOWED_ATOMIC_SYMBOLS:
                    is_valid_composition = False
                    break
            
            if not is_valid_composition:
                continue # Rejette le fragment s'il contient un atome interdit (ex: B, Si, Fe...)

            # 3. Filtre Solvants / Ions purs (Sécurité supplémentaire)
            # Même si la composition est bonne, on vérifie les motifs connus
            temp_smi = Chem.MolToSmiles(frag, isomericSmiles=True) #
            temp_mol = Chem.MolFromSmiles(temp_smi)                #
            smi_frag = Chem.MolToSmiles(temp_mol, isomericSmiles=True) #
            # Glycérol (OCC(O)CO) - 6 atomes, composition C/H/O valide, donc faut le bloquer ici
            if smi_frag == "OCC(O)CO":
                continue
            # DMSO, Ethanol, etc. sont déjà bloqués par MIN_HEAVY_ATOMS=6, mais au cas où :
            if smi_frag in ["CS(=O)C", "CCO", "CO", "CC(C)=O"]:
                continue

            candidates.append((num_heavy, frag))

        if not candidates:
            return None

        # 4. Sélection du plus gros fragment valide
        candidates.sort(key=lambda x: x[0], reverse=True)
        best_frag = candidates[0][1]

        # 5. Canonicalisation
        temp_smi = Chem.MolToSmiles(best_frag, isomericSmiles=True) #
        temp_mol = Chem.MolFromSmiles(temp_smi)                     #
        final_smi = Chem.MolToSmiles(temp_mol, isomericSmiles=True) #
        # Double vérification finale (au cas où la canonicalisation aurait changé quelque chose)
        # (Rare, mais prudent)
        final_mol_check = Chem.MolFromSmiles(final_smi)
        if not final_mol_check:
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
rejection_reasons = {"too_small": 0, "bad_atoms": 0, "solvent": 0, "invalid": 0}

print(f"🚀 Démarrage du nettoyage STRICT pour DeepCOY dans : {ROOT}")
print(f"🛡️ Atomes autorisés : {', '.join(sorted(ALLOWED_ATOMIC_SYMBOLS))}\n")

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
            
            # Petit hack pour debugger si besoin : on pourrait retourner la raison du rejet
            # Mais ici on fait simple : None = rejeté
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
        print(f"✅ {sub}: {kept} ligands conservés ({local_rejected} rejetés)")
        total_cleaned += kept
        total_rejected += local_rejected
    else:
        if os.path.exists(output_file):
            os.remove(output_file)
        print(f"❌ {sub}: 0 ligand valide (Fichier supprimé)")
        total_rejected += local_rejected

    total_processed += 1

print(f"\n{'='*50}")
print(f"🏁 TERMINÉ")
print(f"Dossiers traités : {total_processed}")
print(f"Ligands propres générés : {total_cleaned}")
print(f"Lignes rejetées (ions/solvants/atomes interdits) : {total_rejected}")
print(f"📂 Tes fichiers sont maintenant 100% compatibles DeepCOY.")
