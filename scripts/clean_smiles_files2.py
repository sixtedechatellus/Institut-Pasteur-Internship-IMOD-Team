from rdkit import Chem
import os

ROOT = "data/ValidationEcoli_sansKiKd"  # ou "data/ValidationEcoli_sansKiKd"

# Symboles des atomes qu'on considère comme des ions/métaux à exclure s'ils sont seuls
ION_SYMBOLS = {'H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 
               'Ra', 'Sc', 'Y', 'La', 'Ac', 'Ti', 'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Cr', 
               'Mo', 'W', 'Mn', 'Tc', 'Re', 'Fe', 'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 
               'Pd', 'Pt', 'Cu', 'Ag', 'Au', 'Zn', 'Cd', 'Hg', 'Al', 'Ga', 'In', 'Tl', 
               'Sn', 'Pb', 'Bi', 'Po', 'At', 'F', 'Cl', 'Br', 'I', 'At', 'He', 'Ne', 
               'Ar', 'Kr', 'Xe', 'Rn', 'B', 'Si', 'Ge', 'As', 'Sb', 'Te', 'Se'}

def clean_and_filter_smiles(line):
    """
    1. Extrait le SMILES brut.
    2. Sépare les fragments (sels).
    3. Garde uniquement les fragments organiques (non-métalliques).
    4. Retourne le SMILES canonique combiné ou None si tout est ionique.
    """
    # Nettoyage basique de la ligne
    raw_smi = line.split("\t")[0].split()[0]
    
    try:
        mol = Chem.MolFromSmiles(raw_smi)
        if not mol:
            return None
        
        # Séparer en fragments déconnectés (ex: Cl.CC(=O)O -> [Cl], [CC(=O)O])
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        
        organic_frags = []
        
        for frag in frags:
            # Vérifier si le fragment est un "ion/métal pur"
            # Un fragment est considéré comme ionique si TOUS ses atomes lourds sont dans la liste ION_SYMBOLS
            atoms = [a.GetSymbol() for a in frag.GetAtoms()]
            if not atoms:
                continue
            
            is_pure_ion = all(symbol in ION_SYMBOLS for symbol in atoms)
            
            # On garde le fragment s'il contient au moins un atome NON ionique (C, N, O, S, P...)
            if not is_pure_ion:
                organic_frags.append(frag)
        
        if not organic_frags:
            # La molécule était composée à 100% d'ions/métaux -> On rejette tout
            return None
        
        # Combiner les fragments organiques restants en un seul SMILES
        # (Ex: si on avait Cl.CC(=O)O, on garde juste CC(=O)O)
        final_mol = organic_frags[0]
        for i in range(1, len(organic_frags)):
            final_mol = Chem.CombineMols(final_mol, organic_frags[i])
            
        return Chem.MolToSmiles(final_mol)

    except Exception as e:
        return None

# --- MAIN LOOP ---
count_total = 0
count_skipped = 0

for sub in sorted(os.listdir(ROOT)):
    subdir = os.path.join(ROOT, sub)
    if not os.path.isdir(subdir):
        continue
        
    f_in = os.path.join(subdir, "actives_final.ism")
    # On écrase le fichier de sortie (ou on peut changer le nom si on veut une backup)
    f_out = f_in 
    
    if not os.path.exists(f_in): 
        continue

    valid_smiles = []
    
    with open(f_in) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            
            cleaned_smi = clean_and_filter_smiles(line)
            
            if cleaned_smi:
                valid_smiles.append(cleaned_smi)
            else:
                count_skipped += 1
                # Optionnel : voir quels ligands sont jetés
                # print(f"⚠️ Ignoré (ionique/invalide) dans {sub}: {line[:50]}...")

    if valid_smiles:
        # Suppression des doublons de SMILES identiques (optionnel mais recommandé)
        unique_smiles = list(set(valid_smiles))
        
        with open(f_out, "w") as fout:
            for s in unique_smiles:
                fout.write(s + "\n")
        
        print(f"✅ {sub}: {len(unique_smiles)} ligands organiques conservés (sur {len(valid_smiles)} bruts)")
        count_total += len(unique_smiles)
    else:
        print(f"❌ {sub}: Aucun ligand organique valide trouvé (tout était des ions ou invalide)")
        # Optionnel : supprimer le fichier vide
        if os.path.exists(f_out):
            os.remove(f_out)

print(f"\n🏁 Terminé. Total ligands : {count_total} | Ignorés (ions) : {count_skipped}")
