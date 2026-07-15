import pandas as pd

# Configuration
NEW_BIOLIP_FILE = "BioLiP.txt"  # Le fichier .txt téléchargé
OLD_CSV_FILE = "pdb_extracted.csv"   # Votre CSV de référence
OUTPUT_FILE = "Biolip_new_vs_Biolp_article.csv"

# Colonnes BIoLiP2 (21 colonnes selon documentation)
COLUMNS = [
    'PDB_ID', 'Receptor_chain', 'Resolution', 'Binding_site_code',
    'Ligand_ID', 'Ligand_chain', 'Ligand_serial',
    'Binding_site_residues_PDB', 'Binding_site_residues_renum',
    'Catalytic_site_PDB', 'Catalytic_site_renum',
    'EC_number', 'GO_terms',
    'Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB',
    'UniProt_ID', 'PubMed_ID', 'Ligand_residue_seq', 'Receptor_sequence'
]

print("Chargement du nouveau fichier BIoLiP2...")
df_new = pd.read_csv(NEW_BIOLIP_FILE, sep='\t', header=None, names=COLUMNS, low_memory=False)

print(f"✅ Nouveau: {len(df_new)} interactions ({df_new['PDB_ID'].nunique()} PDBs uniques)")

print("\nChargement de l'ancien CSV...")
df_old = pd.read_csv(OLD_CSV_FILE)

# Détection auto de la colonne PDB dans l'ancien fichier
pdb_col_old = None
for col in df_old.columns:
    if 'pdb' in col.lower():
        pdb_col_old = col
        break

if not pdb_col_old:
    raise ValueError(f"Colonne PDB non trouvée. Colonnes disponibles: {list(df_old.columns)}")

# Normalisation (majuscules, sans espaces)
df_new['PDB_ID_CLEAN'] = df_new['PDB_ID'].str.upper().str.strip()
df_old['PDB_ID_CLEAN'] = df_old[pdb_col_old].astype(str).str.upper().str.strip()

print(f"✅ Ancien: {len(df_old)} entrées ({df_old['PDB_ID_CLEAN'].nunique()} PDBs uniques)")

# Identification des nouveaux
old_pdbs = set(df_old['PDB_ID_CLEAN'].unique())
new_pdbs = set(df_new['PDB_ID_CLEAN'].unique())

nouveaux_pdbs = new_pdbs - old_pdbs
print(f"\n🆕 Nouveaux PDBs: {len(nouveaux_pdbs)}")

if len(nouveaux_pdbs) == 0:
    print("Aucun nouveau PDB trouvé.")
else:
    # Filtrage et export
    df_nouveaux = df_new[df_new['PDB_ID_CLEAN'].isin(nouveaux_pdbs)].copy()
    df_nouveaux = df_nouveaux.drop(columns=['PDB_ID_CLEAN'])  # Nettoyage colonne temp
    
    #df_nouveaux.to_csv(OUTPUT_FILE, index=False)
    print(f"💾 Exporté: {OUTPUT_FILE}")
    print(f"   {len(df_nouveaux)} interactions pour {len(nouveaux_pdbs)} nouveaux PDBs")
    print(f"\nExemples de nouveaux PDBs: {list(nouveaux_pdbs)[:5]}")

#%% Application de filtres comme l'article

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_affinities(df):
    """Convertit les affinités string en pKi/pKd numérique."""
    
    def extract_paffinity(aff_str):
        if pd.isna(aff_str) or aff_str == '' or str(aff_str).upper() == 'NAN':
            return np.nan
        
        # Cherche valeur + unité
        match = re.search(r'(\d+\.?\d*)\s*(nM|uM|µM|mM|pM|fM)', str(aff_str), re.IGNORECASE)
        if not match:
            return np.nan
            
        val, unit = float(match.group(1)), match.group(2).lower()
        
        # Conversion en M puis pKi (en supposant T=298K, mais on prend -log10 pour comparaison)
        conversions = {'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'µm': 1e-6, 'mm': 1e-3}
        molar = val * conversions.get(unit, 1e-9)
        
        # pKi = -log10(Ki[M])
        return -np.log10(molar) if molar > 0 else np.nan
    
    # Prend la meilleure affinité disponible parmi les 4 sources
    aff_cols = ['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']
    
    for col in aff_cols:
        df[f'{col}_pKi'] = df[col].apply(extract_paffinity)
    
    # Meilleure pKi (la plus élevée = affinité la plus forte)
    df['pAffinity'] = df[[f'{c}_pKi' for c in aff_cols]].max(axis=1, skipna=True)
    
    # Pour le tri : si NaN, mettre -inf (sera éliminé par meilleure affinité)
    df['pAffinity_rank'] = df['pAffinity'].fillna(-999)
    
    print(f"Affinités parsées : {df['pAffinity'].notna().sum()}/{len(df)} entrées avec valeur")
    return df

def cluster_sequences_cdhit(df, sequence_col='Receptor_sequence', threshold=0.9):
    """
    Clusterise les séquences à threshold% d'identité avec CD-HIT.
    Retourne un dict {cluster_id: [liste_des_indices]}
    """
    # Vérifier si cd-hit est dispo
    try:
        result = subprocess.run(['cd-hit', '-h'], capture_output=True, text=True)
        cdhit_available = True
    except FileNotFoundError:
        cdhit_available = False
        print("⚠️  CD-HIT non trouvé. Utilisation du fallback UniPROT ID (exact match uniquement).")
        print("   Pour 90% identity, installez CD-HIT : conda install -c bioconda cd-hit")
    
    # Créer fichier FASTA temporaire
    # ID format : index_original|Uniprot pour traçabilité
    records = []
    valid_indices = []
    
    for idx, row in df.iterrows():
        seq = str(row[sequence_col]).strip()
        if len(seq) > 20:  # Ignorer trop courts
            rec = SeqRecord(
                Seq(seq),
                id=f"idx{idx}|{row.get('UniProt_ID', 'unknown')}",
                description=""
            )
            records.append(rec)
            valid_indices.append(idx)
    
    if len(records) == 0:
        print("Aucune séquence valide trouvée")
        return {i: [i] for i in df.index}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_in = os.path.join(tmpdir, "input.fasta")
        fasta_out = os.path.join(tmpdir, "output")
        
        SeqIO.write(records, fasta_in, "fasta")
        
        if cdhit_available:
            # Paramètres CD-HIT : -c 0.9 (90%), -n 5 (word size pour 90%), -d 0 (descriptif complet)
            cmd = [
                'cd-hit', 
                '-i', fasta_in, 
                '-o', fasta_out, 
                '-c', str(threshold),
                '-n', '5',
                '-d', '0',
                '-s', '0.9',  # longueur similaire
                '-S', '0'     # tolerance len
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Erreur CD-HIT: {result.stderr}")
                return {i: [i] for i in valid_indices}
            
            # Parser le fichier .clstr
            clstr_file = fasta_out + ".clstr"
            return parse_cdhit_clusters(clstr_file, df.index)
        else:
            # Fallback : clustering par UniProt ID (100% identité)
            print("Fallback : regroupement par UniProt_ID...")
            groups = df.groupby('UniProt_ID').apply(lambda x: list(x.index)).to_dict()
            return {f"cluster_{k}": v for k, v in groups.items() if k != 'unknown'}

def parse_cdhit_clusters(clstr_file, all_indices):
    """Parse le fichier de clusters CD-HIT."""
    clusters = {}
    current_cluster = None
    
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = line[1:].replace(' ', '_')
                clusters[current_cluster] = []
            else:
                # Format: 0	708aa, >idx123|P00442... at 90.00%
                match = re.search(r'>idx(\d+)\|', line)
                if match:
                    idx = int(match.group(1))
                    if current_cluster:
                        clusters[current_cluster].append(idx)
    
    # Filtrer les clusters vides et convertir index
    valid_clusters = {k: [i for i in v if i in all_indices] 
                     for k, v in clusters.items() if v}
    
    print(f"Clustering terminé : {len(valid_clusters)} clusters créés")
    return valid_clusters

def filter_protein_redundancy(df, clusters):
    """
    Pour chaque cluster de séquences (>90% identité), 
    garde l'entrée avec la meilleure affinité.
    """
    indices_to_keep = []
    
    for cluster_id, idx_list in clusters.items():
        if len(idx_list) == 1:
            indices_to_keep.extend(idx_list)
        else:
            sub_df = df.loc[idx_list]
            # Max pAffinity (si égalité, prendre meilleure résolution ou premier)
            best_idx = sub_df['pAffinity_rank'].idxmax()
            indices_to_keep.append(best_idx)
    
    result = df.loc[indices_to_keep].copy()
    print(f"Filtrage redondance : {len(df)} → {len(result)} entrées ({len(df)-len(result)} supprimées)")
    return result

#def filter_ligand_promiscuity(df):
    """
    Pour chaque ligand apparaissant avec plusieurs cibles,
    garde uniquement la paire avec la meilleure affinité.
    """
    # Grouper par Ligand_ID
#    grouped = df.groupby('Ligand_ID')
    
    # Pour chaque groupe, garder l'index du max pAffinity_rank
#    best_indices = grouped['pAffinity_rank'].idxmax()
    
#    result = df.loc[best_indices].copy()
    
    # Stats
#    promiscuous = grouped.size()
#    n_promiscuous = (promiscuous > 1).sum()
#    n_removed = len(df) - len(result)
    
#    print(f"Filtrage promiscuité : {n_promiscuous} ligands sur {len(grouped)} étaient promiscus")
#    print(f"                     : {len(df)} → {len(result)} entrées ({n_removed} supprimées)")
    
#    return result

def main_pipeline(input_csv, output_csv):
    print("="*60)
    print("FILTRAGE BIoLiP2 : Redondance + Sélectivité")
    print("="*60)
    
    # 1. Chargement
    print(f"\nChargement de {input_csv}...")
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"Entrées initiales : {len(df)}")
    
    # 2. Parsing affinités
    df = parse_affinities(df)
    
    # 3. Clustering séquences (90% identité)
    print("\nClustering des séquences à 90% d'identité...")
    clusters = cluster_sequences_cdhit(df, threshold=0.9)
    
    # 4. Filtre redondance protéique
    df_final = filter_protein_redundancy(df, clusters)
    
    # 5. Filtre promiscuité ligand
#    df_final = filter_ligand_promiscuity(df_step1)
    
    # 6. Export
    # Nettoyage des colonnes auxiliaires
    cols_to_drop = [c for c in df_final.columns if c.endswith('_pKi') or c in ['pAffinity_rank']]
    df_final = df_final.drop(columns=cols_to_drop, errors='ignore')
    
    df_final.to_csv(output_csv, index=False)
    
    print("\n" + "="*60)
    print(f"RÉSULTAT FINAL : {len(df_final)} paires protéine-ligand")
    print(f"Fichier sauvegardé : {output_csv}")
    print("="*60)
    
    # Résumé qualité
    strong = (df_final['pAffinity'] >= 6).sum()  # pKi >= 6 ( Ki <= 1µM)
    print(f"Dont {strong} avec affinité forte (pKi ≥ 6)")

if __name__ == "__main__":
    # Config
    INPUT = "Résultats/Résultats finaux propres/Biolip_new_vs_Biolp_article.csv"  # Le fichier de l'étape précédente
    OUTPUT = "biolip_filtered_90identity.csv"
    
    main_pipeline(INPUT, OUTPUT)

#%% Observation des affinités

import pandas as pd
import numpy as np

df = pd.read_csv("Résultats/Résultats finaux propres/Biolip2_2026_vs_biolip2_article_plusieurs_ligands.csv")

# Analyse détaillée des affinités
total = len(df)
with_affinity = df['pAffinity'].notna().sum()
without_affinity = df['pAffinity'].isna().sum()

print(f"Total paires : {total}")
print(f"Avec affinité mesurée : {with_affinity} ({with_affinity/total*100:.1f}%)")
print(f"Sans affinité (NaN) : {without_affinity} ({without_affinity/total*100:.1f}%)")

if with_affinity > 0:
    print(f"\nRépartition des affinités mesurées :")
    print(f"  pKi ≥ 8 (très forte) : {(df['pAffinity'] >= 8).sum()}")
    print(f"  6 ≤ pKi < 8 (forte) : {((df['pAffinity'] >= 6) & (df['pAffinity'] < 8)).sum()}")
    print(f"  5 ≤ pKi < 6 (modérée) : {((df['pAffinity'] >= 5) & (df['pAffinity'] < 6)).sum()}")
    print(f"  < 5 (faible) : {(df['pAffinity'] < 5).sum()}")

# Vérification : les 6051 correspondent à quoi exactement ?
strong_binders = (df['pAffinity'] >= 6).sum()
others = total - strong_binders
print(f"\nDonc vos {others} paires restantes comprennent :")
print(f"  - Affinité faible/modérée : {with_affinity - strong_binders}")
print(f"  - Aucune mesure : {without_affinity}")


#%% Vérification cofacteurs


import pandas as pd

df = pd.read_csv("biolip_filtered_90identity.csv")

# Liste de cofacteurs / molécules non-drug à exclure
COFACTORS = {
    'NAD', 'NAP', 'FAD', 'FMN', 'COA', 'SAM', 'SAH', 
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',
    'HEM', 'HEC', 'PLP', 'TPP', 'B12',
    'MG', 'ZN', 'CA', 'FE', 'MN', 'CU',
    'SO4', 'PO4', 'GOL', 'EDO', 'DMS', 'ACT'
}

# Identifier les cofacteurs
if 'Ligand_ID' in df.columns:
    df['is_cofactor'] = df['Ligand_ID'].str.upper().isin(COFACTORS)
    
    n_cofactors = df['is_cofactor'].sum()
    print(f"Cofacteurs / non-drugs détectés : {n_cofactors}/{len(df)}")
    
    if n_cofactors > 0:
        print("\nCofacteurs trouvés :")
        print(df[df['is_cofactor']]['Ligand_ID'].value_counts())
    
    # Garder uniquement les vrais inhibiteurs
    df_clean = df[~df['is_cofactor']].copy()
    print(f"\nAprès nettoyage : {len(df_clean)} paires")
    df_clean.to_csv("Biolip2_2026_vs_biolip2_article_plusieurs_ligands.csv", index=False)

#%% Version des filtres comme l'article avec filtre Ki/Kd avant

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def filter_ki_kd_strict(df, aff_cols=['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']):
    """
    Filtre strict : garde uniquement les lignes avec Ki ou Kd explicites.
    Supprime IC50, EC50, et valeurs non spécifiées AVANT tout traitement.
    """
    print(f"\n🔬 Filtre Ki/Kd strict sur {len(df)} entrées...")
    
    def has_ki_or_kd(row):
        """Vérifie si au moins une colonne contient Ki ou Kd (et pas IC50/EC50)."""
        for col in aff_cols:
            if col not in row or pd.isna(row[col]):
                continue
                
            val = str(row[col]).upper()
            
            # Doit contenir KI ou KD explicitement
            has_ki_kd = bool(re.search(r'\bKI\b|\bKD\b', val))
            
            # Ne doit PAS contenir IC50 ou EC50
            has_ic50_ec50 = bool(re.search(r'\bIC50\b|\bEC50\b', val))
            
            if has_ki_kd and not has_ic50_ec50:
                # Extraction rapide pour vérifier que c'est une valeur numérique valide
                if re.search(r'\d+\.?\d*\s*(pM|nM|µM|uM|mM|fM)', val, re.IGNORECASE):
                    return True
        return False
    
    # Application du masque
    mask = df.apply(has_ki_or_kd, axis=1)
    df_filtered = df[mask].copy()
    
    n_removed = len(df) - len(df_filtered)
    print(f"   ✅ Conservé : {len(df_filtered)} entrées avec Ki/Kd")
    print(f"   ❌ Retiré : {n_removed} entrées (IC50, EC50, ou type non spécifié)")
    
    if len(df_filtered) == 0:
        print("⚠️  ATTENTION : Aucune entrée avec Ki/Kd trouvée !")
        print("   Colonnes disponibles :", df.columns.tolist())
        print("   Exemple valeurs :", df[aff_cols[0]].head(3).tolist() if aff_cols[0] in df.columns else "N/A")
    
    return df_filtered

def parse_affinities(df, aff_cols=['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']):
    """Convertit UNIQUEMENT les affinités Ki/Kd en pKi/pKd."""
    
    def extract_paffinity(aff_str):
        if pd.isna(aff_str) or aff_str == '' or str(aff_str).upper() == 'NAN':
            return np.nan
        
        # Vérification préliminaire : on ne parse que si c'est du Ki/Kd
        val_upper = str(aff_str).upper()
        if not (re.search(r'\bKI\b|\bKD\b', val_upper) and not re.search(r'\bIC50\b|\bEC50\b', val_upper)):
            return np.nan
        
        # Cherche valeur + unité
        match = re.search(r'(\d+\.?\d*[eE]?[+-]?\d*)\s*(pM|nM|uM|µM|mM|fM)', str(aff_str), re.IGNORECASE)
        if not match:
            return np.nan
            
        val, unit = float(match.group(1)), match.group(2).lower()
        
        # Conversion en M puis pKi
        conversions = {'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'µm': 1e-6, 'mm': 1e-3, 'fm': 1e-15}
        molar = val * conversions.get(unit, np.nan)
        
        if pd.isna(molar) or molar <= 0:
            return np.nan
            
        return -np.log10(molar)
    
    # Parse chaque colonne
    valid_cols = [c for c in aff_cols if c in df.columns]
    for col in valid_cols:
        df[f'{col}_pKi'] = df[col].apply(extract_paffinity)
    
    # Meilleure pKi parmi les sources disponibles
    pki_cols = [f'{c}_pKi' for c in valid_cols]
    df['pAffinity'] = df[pki_cols].max(axis=1, skipna=True)
    
    # Pour le tri
    df['pAffinity_rank'] = df['pAffinity'].fillna(-999)
    
    # Stats par type (affichage)
    valid = df['pAffinity'].notna().sum()
    print(f"\n📊 Affinités parsées : {valid}/{len(df)} entrées avec valeur numérique")
    
    # Distribution qualitative
    if valid > 0:
        strong = (df['pAffinity'] >= 6).sum()
        medium = ((df['pAffinity'] >= 5) & (df['pAffinity'] < 6)).sum()
        weak = (df['pAffinity'] < 5).sum()
        print(f"   💎 Forte (pKi≥6) : {strong} | Moyenne (5-6) : {medium} | Faible (<5) : {weak}")
    
    return df

def cluster_sequences_cdhit(df, sequence_col='Receptor_sequence', threshold=0.9):
    """Clusterise les séquences à threshold% d'identité avec CD-HIT."""
    try:
        subprocess.run(['cd-hit', '-h'], capture_output=True, text=True)
        cdhit_available = True
    except FileNotFoundError:
        cdhit_available = False
        print("⚠️  CD-HIT non trouvé. Fallback sur séquence exacte.")
    
    records = []
    valid_indices = []
    
    # Détection auto de la colonne séquence
    if sequence_col not in df.columns:
        seq_candidates = [c for c in df.columns if 'seq' in c.lower()]
        if seq_candidates:
            sequence_col = seq_candidates[0]
    
    for idx, row in df.iterrows():
        seq = str(row.get(sequence_col, '')).strip().replace(' ', '').replace('-', '')
        if len(seq) > 20:
            rec = SeqRecord(Seq(seq), id=f"idx{idx}", description="")
            records.append(rec)
            valid_indices.append(idx)
    
    if len(records) == 0:
        return {i: [i] for i in df.index}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_in = os.path.join(tmpdir, "input.fasta")
        fasta_out = os.path.join(tmpdir, "output")
        
        SeqIO.write(records, fasta_in, "fasta")
        
        if cdhit_available:
            cmd = ['cd-hit', '-i', fasta_in, '-o', fasta_out, 
                   '-c', str(threshold), '-n', '5', '-d', '0', '-M', '0', '-T', '0']
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                return {i: [i] for i in valid_indices}
            
            return parse_cdhit_clusters(fasta_out + ".clstr", valid_indices)
        else:
            # Fallback: regroupement exact par séquence
            groups = df.groupby(sequence_col).apply(lambda x: list(x.index)).to_dict()
            return {f"cluster_{k}": v for k, v in groups.items()}

def parse_cdhit_clusters(clstr_file, valid_indices):
    """Parse le fichier .clstr de CD-HIT."""
    clusters = {}
    current_cluster = None
    
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = line[1:].replace(' ', '_')
                clusters[current_cluster] = []
            else:
                match = re.search(r'>idx(\d+)', line)
                if match and current_cluster:
                    idx = int(match.group(1))
                    if idx in valid_indices:
                        clusters[current_cluster].append(idx)
    
    valid_clusters = {k: v for k, v in clusters.items() if v}
    print(f"✅ Clustering terminé : {len(valid_clusters)} clusters")
    return valid_clusters

def filter_protein_redundancy(df, clusters):
    """Pour chaque cluster, garde l'entrée avec la meilleure affinité."""
    indices_to_keep = []
    
    for cluster_id, idx_list in clusters.items():
        valid_idx = [i for i in idx_list if i in df.index]
        if not valid_idx:
            continue
            
        if len(valid_idx) == 1:
            indices_to_keep.extend(valid_idx)
        else:
            sub_df = df.loc[valid_idx]
            best_idx = sub_df['pAffinity_rank'].idxmax()
            indices_to_keep.append(best_idx)
    
    result = df.loc[indices_to_keep].copy()
    print(f"🧬 Filtre redondance : {len(df)} → {len(result)} ({len(df)-len(result)} doublons)")
    return result

#def filter_ligand_promiscuity(df, ligand_col='Ligand_ID'):
    """Pour chaque ligand, garde la paire avec la meilleure affinité."""
#    if ligand_col not in df.columns:
#        print(f"⚠️ Colonne {ligand_col} non trouvée")
#        return df
    
#    grouped = df.groupby(ligand_col)
#    best_indices = grouped['pAffinity_rank'].idxmax()
    
#    result = df.loc[best_indices].copy()
#    n_promiscuous = (grouped.size() > 1).sum()
    
#    print(f"🧪 Filtre promiscuité : {n_promiscuous} ligands promiscus traités")
#    print(f"                    : {len(df)} → {len(result)} ({len(df)-len(result)} retirés)")
    
#    return result

def main_pipeline(input_csv, output_csv):
    print("="*60)
    print("FILTRAGE BIoLiP2 : Ki/Kd Strict → Redondance → Sélectivité")
    print("="*60)
    
    # 1. Chargement
    print(f"\n📂 Chargement : {input_csv}")
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"   Entrées initiales : {len(df)}")
    
    # 2. FILTRE STRICT Ki/Kd (NOUVEAU - AVANT tout le reste)
    df = filter_ki_kd_strict(df)
    
    if len(df) == 0:
        print("❌ Arrêt : aucune donnée Ki/Kd exploitable")
        return
    
    # 3. Parsing affinités (maintenant que c'est du Ki/Kd pur)
    df = parse_affinities(df)
    
    # Vérifie qu'on a bien des valeurs numériques
    valid_aff = df['pAffinity'].notna().sum()
    if valid_aff == 0:
        print("❌ Erreur : aucune affinité Ki/Kd parsable numériquement")
        return
    
    # 4. Clustering
    print(f"\n🔄 Clustering séquences...")
    clusters = cluster_sequences_cdhit(df, threshold=0.9)
    
    # 5. Filtre redondance protéique
    df_final = filter_protein_redundancy(df, clusters)
    
    # 6. Filtre promiscuité ligand
#    df_final = filter_ligand_promiscuity(df_step1)
    
    # 7. Export
    cols_to_drop = [c for c in df_final.columns if c.endswith('_pKi') or c == 'pAffinity_rank']
    df_export = df_final.drop(columns=cols_to_drop, errors='ignore')
    
    df_export.to_csv(output_csv, index=False)
    
    # Résumé final
    print("\n" + "="*60)
    print(f"🎯 RÉSULTAT FINAL")
    print(f"   Paires protéine-ligand : {len(df_final)}")
    print(f"   Réduction totale : {len(df)} → {len(df_final)} ({(1-len(df_final)/len(df))*100:.1f}%)")
    print(f"   Qualité : {(df_final['pAffinity'] >= 6).sum()} fortes (pKi≥6)")
    print(f"   Fichier : {output_csv}")
    print("="*60)
    
    # Validation scientifique
    print("\n✅ Validation : Données comparables (toutes en pKi/pKd thermodynamique)")

if __name__ == "__main__":
    INPUT = "Résultats/Résultats finaux propres/Biolip_new_vs_Biolp_article.csv"
    OUTPUT = "biolip_filtered_90identity.csv"
    
    main_pipeline(INPUT, OUTPUT)

#%% FIltre anti leakage (30% Idt avec bolip2 article)

#!/usr/bin/env python3
"""
Vérifie qu'il n'y a pas de fuite (leakage) entre :
  - Dataset de validation (avec UniProt_ID)
  - Dataset de finetuning (PDB IDs → SIFTS → UniProt)

Checks :
  1) PDB identiques
  2) UniProt identiques
  3) Identité de séquence ≥30% (k-mer proxy)
"""

import pandas as pd
import os
import urllib.request
import gzip
import shutil
import time
import requests

# =============================================================
# CONFIG
# =============================================================
VALIDATION_FILE = "Résultats/Résultats finaux propres/BioLiP_ecoli_ki_kd_ic50.csv"       # ← ton dataset validation
FINETUNE_FILE = "PDB finetuning article.csv"     # ← tes PDB de finetuning

# Nom de la colonne PDB dans le fichier finetuning
FINETUNE_PDB_COL = "PDB_ID"  # ou "pdb_ids" — adapter

OUTPUT_DIR = "leakage_check_E.coli"
IDENTITY_THRESHOLD = 0.3

SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
SIFTS_GZ = "pdb_chain_uniprot.tsv.gz"
SIFTS_TSV = "pdb_chain_uniprot.tsv"


# =============================================================
# SIFTS
# =============================================================
def get_sifts():
    if not os.path.exists(SIFTS_TSV):
        if not os.path.exists(SIFTS_GZ):
            print("[SIFTS] Téléchargement (~80 Mo)...")
            urllib.request.urlretrieve(SIFTS_URL, SIFTS_GZ)
        print("[SIFTS] Décompression...")
        with gzip.open(SIFTS_GZ, "rb") as f_in, open(SIFTS_TSV, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        print("[SIFTS] Déjà présent.")

    print("[SIFTS] Chargement...")
    header_line = 0
    with open(SIFTS_TSV) as f:
        for i, line in enumerate(f):
            if line.startswith("PDB"):
                header_line = i
                break
    df = pd.read_csv(SIFTS_TSV, sep="\t", skiprows=header_line)
    df.columns = df.columns.str.strip()
    df["PDB"] = df["PDB"].str.strip().str.lower()
    df["SP_PRIMARY"] = df["SP_PRIMARY"].str.strip()
    print(f"    {len(df)} mappings")
    return df[["PDB", "SP_PRIMARY"]].drop_duplicates()


# =============================================================
# Charger datasets
# =============================================================
def load_validation(filepath):
    print(f"\n[VALIDATION] Chargement : {filepath}")
    df = pd.read_csv(filepath)
    df["_pdb"] = df["PDB_ID"].astype(str).str.strip().str.lower()
    df["_uniprot"] = df["UniProt_ID"].astype(str).str.strip()

    pdbs = set(df["_pdb"].dropna().unique()) - {"nan", ""}
    uniprots = set(df["_uniprot"].dropna().unique()) - {"nan", ""}

    print(f"    {len(df)} lignes")
    print(f"    {len(pdbs)} PDB uniques")
    print(f"    {len(uniprots)} UniProt uniques")
    return df, pdbs, uniprots


def load_finetune(filepath, sifts):
    print(f"\n[FINETUNE] Chargement : {filepath}")
    df = pd.read_csv(filepath)
    df["_pdb"] = df[FINETUNE_PDB_COL].astype(str).str.strip().str.lower()

    pdbs = set(df["_pdb"].dropna().unique()) - {"nan", ""}
    print(f"    {len(df)} PDB")

    # Mapper vers UniProt via SIFTS
    merged = df.merge(sifts, left_on="_pdb", right_on="PDB", how="left")
    uniprots = set(merged["SP_PRIMARY"].dropna().unique())
    print(f"    → {len(uniprots)} UniProt uniques")
    return df, pdbs, uniprots


# =============================================================
# Check 1 & 2 : PDB et UniProt identiques
# =============================================================
def check_exact_overlap(val_pdbs, val_uniprots, ft_pdbs, ft_uniprots):
    print("\n" + "=" * 60)
    print("CHECK 1 : PDB identiques")
    print("=" * 60)

    pdb_overlap = val_pdbs & ft_pdbs
    print(f"  Validation PDB  : {len(val_pdbs)}")
    print(f"  Finetuning PDB  : {len(ft_pdbs)}")
    print(f"  PDB en commun   : {len(pdb_overlap)}")

    if pdb_overlap:
        print(f"\n  ⚠ PDB PARTAGÉS :")
        for p in sorted(list(pdb_overlap)[:20]):
            print(f"    {p}")
        if len(pdb_overlap) > 20:
            print(f"    ... et {len(pdb_overlap) - 20} de plus")
    else:
        print(f"  ✅ Aucun PDB partagé")

    print("\n" + "=" * 60)
    print("CHECK 2 : UniProt identiques")
    print("=" * 60)

    uni_overlap = val_uniprots & ft_uniprots
    print(f"  Validation UniProt  : {len(val_uniprots)}")
    print(f"  Finetuning UniProt  : {len(ft_uniprots)}")
    print(f"  UniProt en commun   : {len(uni_overlap)}")

    if uni_overlap:
        print(f"\n  ⚠ UniProt PARTAGÉS ({len(uni_overlap)}) :")
        for u in sorted(list(uni_overlap)[:30]):
            print(f"    {u}")
        if len(uni_overlap) > 30:
            print(f"    ... et {len(uni_overlap) - 30} de plus")
    else:
        print(f"  ✅ Aucun UniProt partagé")

    return pdb_overlap, uni_overlap


# =============================================================
# Télécharger séquences
# =============================================================
def download_sequences(uniprot_ids, output_fasta, label):
    print(f"\n[DL] Séquences {label} ({len(uniprot_ids)})...")

    if os.path.exists(output_fasta):
        n = sum(1 for l in open(output_fasta) if l.startswith(">"))
        if n >= len(uniprot_ids) * 0.8:
            print(f"    Déjà {n} séquences, skip.")
            return
        else:
            print(f"    Seulement {n}, on recommence...")
            os.remove(output_fasta)

    ids = sorted(uniprot_ids)
    batch_size = 50
    n_ok = 0
    n_fail = 0

    with open(output_fasta, "w") as fout:
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            accessions = ",".join(batch)

            try:
                r = requests.get(
                    "https://rest.uniprot.org/uniprotkb/accessions",
                    params={"accessions": accessions, "format": "fasta"},
                    timeout=60
                )
                if r.ok and ">" in r.text:
                    fout.write(r.text.strip() + "\n")
                    n_ok += r.text.count(">")
                else:
                    # Fallback une par une
                    for acc in batch:
                        try:
                            r2 = requests.get(
                                f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                                timeout=30
                            )
                            if r2.ok and ">" in r2.text:
                                fout.write(r2.text.strip() + "\n")
                                n_ok += 1
                            else:
                                n_fail += 1
                        except:
                            n_fail += 1
                        time.sleep(0.1)
            except Exception as e:
                print(f"    Erreur batch {i}: {e}")
                for acc in batch:
                    try:
                        r2 = requests.get(
                            f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                            timeout=30
                        )
                        if r2.ok and ">" in r2.text:
                            fout.write(r2.text.strip() + "\n")
                            n_ok += 1
                        else:
                            n_fail += 1
                    except:
                        n_fail += 1
                    time.sleep(0.1)

            done = min(i + batch_size, len(ids))
            if done % 500 == 0 or done == len(ids):
                print(f"    {done}/{len(ids)} (ok: {n_ok}, fail: {n_fail})")
            time.sleep(0.2)

    print(f"    → {n_ok} séquences → {output_fasta}")


# =============================================================
# Parser FASTA
# =============================================================
def parse_fasta(filepath):
    seqs = {}
    current_id = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                header = line[1:]
                parts = header.split("|")
                current_id = parts[1] if len(parts) >= 2 else header.split()[0]
                current_seq = []
            elif current_id:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs


# =============================================================
# Check 3 : Identité séquence ≥30% par k-mer
# =============================================================
def check_sequence_identity(val_seqs, ft_seqs, threshold=0.3):
    print("\n" + "=" * 60)
    print(f"CHECK 3 : Identité séquence ≥{threshold*100:.0f}%")
    print("=" * 60)

    K = 3
    kmer_threshold = threshold * 0.7  # conservateur

    # Pre-compute k-mers finetuning
    print(f"  Indexation finetuning ({len(ft_seqs)} séquences)...")
    ft_kmers = {}
    for acc, seq in ft_seqs.items():
        if len(seq) >= K:
            ft_kmers[acc] = set(seq[i:i+K] for i in range(len(seq) - K + 1))

    print(f"  Comparaison {len(val_seqs)} validation vs {len(ft_kmers)} finetuning...")

    matches = {}  # val_acc → (best_ft_acc, best_score)
    t0 = time.time()

    for idx, (val_acc, val_seq) in enumerate(val_seqs.items()):
        if len(val_seq) < K:
            continue

        val_kset = set(val_seq[i:i+K] for i in range(len(val_seq) - K + 1))
        n_val = len(val_kset)

        best_score = 0
        best_target = None

        for ft_acc, ft_kset in ft_kmers.items():
            shared = len(val_kset & ft_kset)
            denom = min(n_val, len(ft_kset))
            if denom > 0:
                score = shared / denom
                if score > best_score:
                    best_score = score
                    best_target = ft_acc

        if best_score >= kmer_threshold:
            matches[val_acc] = (best_target, best_score)

        if (idx + 1) % 100 == 0 or idx + 1 == len(val_seqs):
            elapsed = time.time() - t0
            print(f"    {idx+1}/{len(val_seqs)}  leaks={len(matches)}  ({elapsed:.0f}s)")

    # Résultats
    n_clean = len(val_seqs) - len(matches)
    print(f"\n  Résultats :")
    print(f"    Validation séquences         : {len(val_seqs)}")
    print(f"    Avec homologue ≥30% dans FT  : {len(matches)} ⚠ LEAKAGE")
    print(f"    Sans homologue (clean)       : {n_clean} ✅")

    if matches:
        print(f"\n  Détail des leaks (top 30) :")
        sorted_matches = sorted(matches.items(), key=lambda x: -x[1][1])
        for val_acc, (ft_acc, score) in sorted_matches[:30]:
            print(f"    {val_acc} ↔ {ft_acc}  (k-mer score: {score:.2f})")

    return matches


# =============================================================
# Export
# =============================================================
def export_results(pdb_overlap, uni_overlap, seq_matches, val_df):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1) Résumé
    p1 = os.path.join(OUTPUT_DIR, "leakage_summary.txt")
    with open(p1, "w") as f:
        f.write(f"=== LEAKAGE CHECK ===\n\n")
        f.write(f"PDB partagés           : {len(pdb_overlap)}\n")
        f.write(f"UniProt partagés       : {len(uni_overlap)}\n")
        f.write(f"Homologues ≥30% seq    : {len(seq_matches)}\n")
        f.write(f"\nTotal problématiques   : {len(uni_overlap | set(seq_matches.keys()))}\n")

    # 2) Liste des protéines à problème
    all_leaky = uni_overlap | set(seq_matches.keys())
    p2 = os.path.join(OUTPUT_DIR, "leaky_uniprots.txt")
    with open(p2, "w") as f:
        f.write("uniprot_id\tleak_type\tbest_match\tscore\n")
        for u in sorted(all_leaky):
            if u in uni_overlap and u in seq_matches:
                ltype = "exact+sequence"
                mt, sc = seq_matches[u]
                f.write(f"{u}\t{ltype}\t{mt}\t{sc:.3f}\n")
            elif u in uni_overlap:
                f.write(f"{u}\texact_uniprot\t-\t1.0\n")
            else:
                mt, sc = seq_matches[u]
                f.write(f"{u}\tsequence_homolog\t{mt}\t{sc:.3f}\n")

    # 3) Dataset validation nettoyé (sans les leaky)
    val_clean = val_df[~val_df["_uniprot"].isin(all_leaky)].copy()
    p3 = os.path.join(OUTPUT_DIR, "validation_clean.csv")
    # Drop colonnes internes
    cols_to_drop = [c for c in val_clean.columns if c.startswith("_")]
    val_clean.drop(columns=cols_to_drop, inplace=True)
    val_clean.to_csv(p3, index=False)

    # 4) Dataset validation problématique
    val_leaky = val_df[val_df["_uniprot"].isin(all_leaky)].copy()
    p4 = os.path.join(OUTPUT_DIR, "validation_leaky.csv")
    cols_to_drop = [c for c in val_leaky.columns if c.startswith("_")]
    val_leaky.drop(columns=cols_to_drop, inplace=True)
    val_leaky.to_csv(p4, index=False)

    print(f"\n[EXPORT]")
    print(f"  {p1}")
    print(f"  {p2}  ({len(all_leaky)} protéines à problème)")
    print(f"  {p3}  ({len(val_clean)} entrées propres)")
    print(f"  {p4}  ({len(val_leaky)} entrées à retirer)")


# =============================================================
# MAIN
# =============================================================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    val_fasta = os.path.join(OUTPUT_DIR, "validation_seqs.fasta")
    ft_fasta = os.path.join(OUTPUT_DIR, "finetune_seqs.fasta")

    # 1) SIFTS
    sifts = get_sifts()

    # 2) Charger
    val_df, val_pdbs, val_uniprots = load_validation(VALIDATION_FILE)
    ft_df, ft_pdbs, ft_uniprots = load_finetune(FINETUNE_FILE, sifts)

    # 3) Check exact
    pdb_overlap, uni_overlap = check_exact_overlap(
        val_pdbs, val_uniprots, ft_pdbs, ft_uniprots
    )

    # 4) Check séquence — seulement pour ceux PAS déjà flaggés
    val_to_check = val_uniprots - uni_overlap  # pas besoin de re-checker les exacts
    ft_to_check = ft_uniprots

    print(f"\n  Séquences à comparer : {len(val_to_check)} validation vs {len(ft_to_check)} finetuning")

    download_sequences(val_to_check, val_fasta, "validation")
    download_sequences(ft_to_check, ft_fasta, "finetuning")

    val_seqs = parse_fasta(val_fasta)
    ft_seqs = parse_fasta(ft_fasta)
    print(f"  Parsées : validation={len(val_seqs)}, finetuning={len(ft_seqs)}")

    if len(val_seqs) == 0 or len(ft_seqs) == 0:
        print("  ⚠ Pas assez de séquences téléchargées !")
        print("  Le check séquence sera incomplet.")
        seq_matches = {}
    else:
        seq_matches = check_sequence_identity(val_seqs, ft_seqs, IDENTITY_THRESHOLD)

    # 5) Résumé final
    all_leaky = uni_overlap | set(seq_matches.keys())
    n_clean = len(val_uniprots) - len(all_leaky)

    print("\n" + "=" * 60)
    print("RÉSUMÉ FINAL")
    print("=" * 60)
    print(f"  Validation total        : {len(val_uniprots)} protéines")
    print(f"  PDB identiques          : {len(pdb_overlap)}")
    print(f"  UniProt identiques      : {len(uni_overlap)}")
    print(f"  Homologues ≥30%         : {len(seq_matches)}")
    print(f"  Total à retirer         : {len(all_leaky)}")
    print(f"  Restant propre          : {n_clean} ✅")
    print(f"  % leakage               : {100*len(all_leaky)/max(len(val_uniprots),1):.1f}%")

    # 6) Export
    export_results(pdb_overlap, uni_overlap, seq_matches, val_df)

    print(f"\nFichiers dans ./{OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
    
    
#%% Visualisation données ki et IC50

import pandas as pd
import numpy as np

# =============================================================
# 1. CHARGEMENT DES DONNÉES
# =============================================================
# Option A : depuis un fichier CSV
df = pd.read_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv")

# Option B : depuis un fichier Excel
# df = pd.read_excel("votre_fichier.xlsx")

# Option C : données d'exemple pour tester
# df = pd.DataFrame({
#     'Molecule':  ['MolA','MolB','MolC','MolD','MolE','MolF','MolG','MolH'],
#     'pKi':       [9.2,   7.5,   5.3,   6.8,   8.1,   4.2,   5.9,   8.5],
#     'IC50_nM':   [0.6,   32,    5000,  158,   8,     63000, 1260,  3.2]
# })

# Afficher un aperçu
print("=" * 60)
print("APERÇU DES DONNÉES")
print("=" * 60)
print(df.head(10))
print(f"\nNombre total de composés : {len(df)}")

# =============================================================
# 2. ANALYSE pKi
# =============================================================
# --- Adaptez le nom de la colonne pKi si nécessaire ---
PKI_COL = 'pKi_Kd'  # <-- Modifier ici si votre colonne s'appelle autrement

print("\n" + "=" * 60)
print("ANALYSE pKi")
print("=" * 60)

seuils_pki = [5, 6, 7, 8, 9]

for seuil in seuils_pki:
    mask = df[PKI_COL] > seuil
    count = mask.sum()
    pct = (count / len(df)) * 100
    print(f"  pKi > {seuil} : {count:>6} composés  ({pct:>6.1f}%)")

# Détail par tranche
print("\n--- Répartition par tranche de pKi ---")
bins_pki = [0, 5, 6, 7, 8, 9, float('inf')]
labels_pki = ['≤5', '5-6', '6-7', '7-8', '8-9', '>9']
df['pKi_tranche'] = pd.cut(df[PKI_COL], bins=bins_pki, labels=labels_pki, right=True)
print(df['pKi_tranche'].value_counts().sort_index())

# =============================================================
# 3. ANALYSE IC50
# =============================================================
# --- Adaptez le nom de la colonne IC50 et l'unité ---
IC50_COL = 'Affinity_BindingDB'  # <-- Modifier ici
UNITE = 'nM'          # 'nM' ou 'uM'

print("\n" + "=" * 60)
print(f"ANALYSE IC50 ({UNITE})")
print("=" * 60)

# --- Seuils IC50 en nM (adapter si en µM) ---
if UNITE == 'nM':
    seuils_ic50 = [1, 10, 100, 1000, 10000]
    # IC50 < seuil  = composé actif
    for seuil in seuils_ic50:
        mask = df[IC50_COL] < seuil
        count = mask.sum()
        pct = (count / len(df)) * 100
        label = f"{seuil} nM"
        if seuil >= 1000:
            label += f" ({seuil/1000:.0f} µM)"
        print(f"  IC50 < {label:>18} : {count:>6} composés  ({pct:>6.1f}%)")

elif UNITE == 'uM':
    seuils_ic50 = [0.001, 0.01, 0.1, 1, 10]
    for seuil in seuils_ic50:
        mask = df[IC50_COL] < seuil
        count = mask.sum()
        pct = (count / len(df)) * 100
        print(f"  IC50 < {seuil:>8} µM : {count:>6} composés  ({pct:>6.1f}%)")

# --- Conversion IC50 → pIC50 ---
print("\n--- Conversion IC50 → pIC50 ---")
if UNITE == 'nM':
    df['pIC50'] = -np.log10(df[IC50_COL] * 1e-9)
elif UNITE == 'uM':
    df['pIC50'] = -np.log10(df[IC50_COL] * 1e-6)

# Analyse pIC50 avec les mêmes seuils que pKi
print("\n--- Analyse pIC50 (mêmes seuils que pKi) ---")
for seuil in seuils_pki:
    mask = df['pIC50'] > seuil
    count = mask.sum()
    pct = (count / len(df)) * 100
    print(f"  pIC50 > {seuil} : {count:>6} composés  ({pct:>6.1f}%)")

# Répartition par tranche de pIC50
df['pIC50_tranche'] = pd.cut(df['pIC50'], bins=bins_pki, labels=labels_pki, right=True)
print("\n--- Répartition par tranche de pIC50 ---")
print(df['pIC50_tranche'].value_counts().sort_index())

# =============================================================
# 4. TABLEAU RÉCAPITULATIF
# =============================================================
print("\n" + "=" * 60)
print("TABLEAU RÉCAPITULATIF")
print("=" * 60)

recap = []
for seuil in seuils_pki:
    row = {
        'Seuil': f'> {seuil}',
        f'Nb pKi > {seuil}': (df[PKI_COL] > seuil).sum(),
        f'Nb pIC50 > {seuil}': (df['pIC50'] > seuil).sum(),
    }
    recap.append(row)

recap_df = pd.DataFrame(recap)
print(recap_df.to_string(index=False))

# =============================================================
# 5. LISTER LES COMPOSÉS ACTIFS (pKi > 8)
# =============================================================
print("\n" + "=" * 60)
print("COMPOSÉS AVEC pKi > 8 (haute affinité)")
print("=" * 60)
actifs = df[df[PKI_COL] > 8].sort_values(PKI_COL, ascending=False)
print(actifs.to_string(index=False))

# =============================================================
# 6. EXPORT DES RÉSULTATS
# =============================================================
df.to_csv("resultats_analyse.csv", index=False)
print("\n✅ Résultats exportés dans 'resultats_analyse.csv'")

#%% Idem

import pandas as pd
import numpy as np

# =============================================================
# 1. CHARGEMENT
# =============================================================
df = pd.read_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv")

# Vérifier le nom exact de la colonne
print("Colonnes disponibles :")
print(df.columns.tolist())

# =============================================================
# 2. ADAPTEZ LE NOM EXACT DE VOTRE COLONNE ICI
# =============================================================
PIC50_COL = 'pIC50'  # <-- modifiez si besoin (ex: 'pIC50', 'pic50', 'pChEMBL Value'...)

# Vérification
print(f"\nType : {df[PIC50_COL].dtype}")
print(f"Aperçu :")
print(df[PIC50_COL].head(10))

# =============================================================
# 3. CONVERTIR EN NUMÉRIQUE (sécurité)
# =============================================================
df[PIC50_COL] = pd.to_numeric(df[PIC50_COL], errors='coerce')
df_clean = df.dropna(subset=[PIC50_COL]).copy()

print(f"\nValeurs valides : {len(df_clean)} / {len(df)}")

# =============================================================
# 4. COMPTAGE PAR SEUILS
# =============================================================
print("\n" + "=" * 60)
print("ANALYSE pIC50")
print("=" * 60)

seuils = [5, 6, 7, 8, 9]

for seuil in seuils:
    count = (df_clean[PIC50_COL] > seuil).sum()
    total = len(df_clean)
    pct = (count / total * 100)
    print(f"  pIC50 > {seuil}  :  {count:>6} / {total}  ({pct:>6.1f}%)")

# =============================================================
# 5. RÉPARTITION PAR TRANCHE
# =============================================================
print("\n--- Répartition par tranche ---")
bins   = [0, 5, 6, 7, 8, 9, float('inf')]
labels = ['≤5', '5-6', '6-7', '7-8', '8-9', '>9']
df_clean['tranche'] = pd.cut(df_clean[PIC50_COL], bins=bins, labels=labels, right=True)
print(df_clean['tranche'].value_counts().sort_index())

# =============================================================
# 6. STATISTIQUES
# =============================================================
print(f"\n--- Statistiques ---")
print(f"  Min    : {df_clean[PIC50_COL].min():.2f}")
print(f"  Max    : {df_clean[PIC50_COL].max():.2f}")
print(f"  Moyenne: {df_clean[PIC50_COL].mean():.2f}")
print(f"  Médiane: {df_clean[PIC50_COL].median():.2f}")

# =============================================================
# 7. TOP 20 COMPOSÉS LES PLUS ACTIFS
# =============================================================
print("\n" + "=" * 60)
print("TOP 20 - MEILLEURS pIC50")
print("=" * 60)
top20 = df_clean.nlargest(20, PIC50_COL)
print(top20.to_string(index=False))