import requests
import pandas as pd
from tqdm import tqdm
from io import StringIO


def get_ecoli_proteome_complete():
    """Récupère le protéome complet d'E. coli K-12 MG1655 depuis UniProt."""
    
    url = "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000625&format=tsv"
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    
    content = bytearray()
    with tqdm(total=total_size, unit='B', unit_scale=True, desc="Download") as pbar:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                content.extend(chunk)
                pbar.update(len(chunk))
    
    df = pd.read_csv(StringIO(content.decode('utf-8')), sep='\t', low_memory=False)
    
    column_mapping = {
        'Entry': 'uniprot_id',
        'Entry Name': 'entry_name',
        'Gene Names': 'gene_names',
        'Protein names': 'protein_name',
        'Organism': 'organism',
        'Length': 'length',
        'Sequence': 'sequence',
        'Cross-reference (PDB)': 'pdb_ids',
        'Cross-reference (AlphaFoldDB)': 'alphafold_ids',
        'Cross-reference (STRING)': 'string_ids',
        'Reviewed': 'reviewed'
    }
    
    existing_cols = {k: v for k, v in column_mapping.items() if k in df.columns}
    df = df.rename(columns=existing_cols)
    
    if 'pdb_ids' in df.columns:
        df['pdb_ids'] = df['pdb_ids'].fillna('')
        df['pdb_ids'] = df['pdb_ids'].str.replace(r'\s*;\s*', ';', regex=True).str.strip(';')
        df['has_pdb'] = df['pdb_ids'].apply(lambda x: len(str(x)) > 0)
    else:
        df['pdb_ids'] = ''
        df['has_pdb'] = False
    
    if 'alphafold_ids' in df.columns:
        df['alphafold_ids'] = df['alphafold_ids'].fillna('')
        df['alphafold_ids'] = df['alphafold_ids'].str.replace(r'\s*;\s*', ';', regex=True).str.strip(';')
        df['has_alphafold'] = df['alphafold_ids'].apply(lambda x: len(str(x)) > 0)
    else:
        df['alphafold_ids'] = ''
        df['has_alphafold'] = False
    
    return df


df = get_ecoli_proteome_complete()
df.to_csv("ecoli_k12_COMPLETE_proteome.csv", index=False)


#%% Query PDB et structures alphafold


import time


def map_structures_final(uniprot_csv, output_csv="ecoli_structures_final.csv"):
    """Mappe les structures PDB et AlphaFold pour chaque protéine E. coli."""
    
    df = pd.read_csv(uniprot_csv)
    uniprot_ids = df['uniprot_id'].tolist()
    
    pdb_results = {}
    alphafold_results = {}
    
    for i in tqdm(range(0, len(uniprot_ids), 100)):
        batch = uniprot_ids[i:i+100]
        ids_str = ",".join(batch)
        url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={ids_str}&fields=accession,xref_pdb,xref_alphafolddb&format=json"
        
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 200:
                data = r.json()
                for entry in data.get("results", []):
                    accession = entry.get("primaryAccession")
                    pdb_entries = []
                    
                    for ref in entry.get("uniProtKBCrossReferences", []):
                        if ref.get("database") == "PDB":
                            pdb_entries.append(ref.get("id"))
                        elif ref.get("database") == "AlphaFoldDB":
                            alphafold_results[accession] = ref.get("id")
                    
                    if pdb_entries:
                        pdb_results[accession] = ";".join(pdb_entries)
                    
                    if accession not in alphafold_results:
                        alphafold_results[accession] = f"AF-{accession}-F1"
            else:
                print(f"Status {r.status_code} for batch {i}")
        except Exception as e:
            print(f"Error batch {i}: {e}")
        
        time.sleep(0.2)
    
    df['pdb_ids'] = df['uniprot_id'].map(pdb_results).fillna('')
    df['has_pdb'] = df['pdb_ids'] != ''
    df['alphafold_ids'] = df['uniprot_id'].map(alphafold_results)
    
    df['structure_id'] = df.apply(
        lambda row: row['pdb_ids'].split(';')[0] if row['has_pdb'] else row['alphafold_ids'],
        axis=1
    )
    df['structure_source'] = df['has_pdb'].apply(lambda x: 'PDB' if x else 'AlphaFold')
    
    df['alphafold_fallback'] = df.apply(
        lambda row: f"AF-{row['uniprot_id']}-F1" if not row['has_pdb'] else '',
        axis=1
    )
    
    df.to_csv(output_csv, index=False)
    return df


df_result = map_structures_final("ecoli_k12_COMPLETE_proteome.csv")


#%% Etude pour voir quelles protéines d'E coli sont dans le dataset de finetuning


def match_with_biolip2(ecoli_csv, biolip2_csv, output_csv="ecoli_in_biolip2.csv"):
    """Matche les PDB du protéome E. coli avec ceux présents dans BioLip2."""
    
    df_ecoli = pd.read_csv(ecoli_csv)
    df_biolip = pd.read_csv(biolip2_csv)
    
    biolip_pdb_col = df_biolip.columns[0]
    biolip_pdbs = set(df_biolip[biolip_pdb_col].astype(str).str.upper().str.strip())
    
    def find_matching_pdbs(pdb_cell):
        if pd.isna(pdb_cell) or pdb_cell == '':
            return []
        pdbs = [p.strip().upper() for p in str(pdb_cell).split(';') if len(p.strip()) == 4]
        return [p for p in pdbs if p in biolip_pdbs]
    
    df_ecoli['matching_pdbs'] = df_ecoli['pdb_ids'].apply(find_matching_pdbs)
    df_ecoli['n_matches'] = df_ecoli['matching_pdbs'].apply(len)
    df_ecoli['has_biolip_match'] = df_ecoli['n_matches'] > 0
    
    df_matches = df_ecoli[df_ecoli['has_biolip_match']].copy()
    df_matches['pdb_ids_original'] = df_matches['pdb_ids']
    df_matches['pdb_ids'] = df_matches['matching_pdbs'].apply(lambda x: ';'.join(x))
    df_matches_clean = df_matches.drop(columns=['matching_pdbs', 'n_matches', 'has_biolip_match'])
    
    df_matches_clean.to_csv(output_csv, index=False)
    return df_matches_clean


if __name__ == "__main__":
    ECOLI_FILE = "Résultats/Résultats finaux propres/Escherichia_Coli_K12_protéome_complet.csv"
    BIOLIP_FILE = "PDB finetuning article.csv"
    OUTPUT = "ecoli_biolip2_intersection.csv"
    
    df_result = match_with_biolip2(ECOLI_FILE, BIOLIP_FILE, OUTPUT)


#%% Protéines avec un PDB mais pas das biolip2 article


def extract_not_in_biolip2(ecoli_csv, biolip2_csv, output_csv="ecoli_pdb_not_in_biolip2.csv"):
    """Extrait les protéines avec PDB expérimental mais absentes de BioLip2."""
    
    df_ecoli = pd.read_csv(ecoli_csv)
    df_biolip = pd.read_csv(biolip2_csv)
    biolip_pdbs = set(df_biolip.iloc[:, 0].astype(str).str.upper().str.strip())
    
    def has_any_match(pdb_cell):
        if pd.isna(pdb_cell) or pdb_cell == '':
            return False
        pdbs = [p.strip().upper() for p in str(pdb_cell).split(';')]
        return any(p in biolip_pdbs for p in pdbs)
    
    df_ecoli['in_biolip2'] = df_ecoli['pdb_ids'].apply(has_any_match)
    mask = (df_ecoli['has_pdb'] == True) & (~df_ecoli['in_biolip2'])
    df_not_in = df_ecoli[mask].copy()
    df_not_in = df_not_in.drop(columns=['in_biolip2'], errors='ignore')
    df_not_in.to_csv(output_csv, index=False)
    
    return df_not_in


df_hors_biolip = extract_not_in_biolip2(
    "Résultats/Résultats finaux propres/Escherichia_Coli_K12_protéome_complet.csv",
    "PDB finetuning article.csv",
    "ecoli_pdb_hors_biolip2.csv"
)


#%% Comparaison avec Biolip2 2026 pour avoir les matchs


def match_ecoli_with_new_biolip2(ecoli_csv, new_biolip_file, output_csv="ecoli_in_new_biolip2.csv"):
    """Matche E. coli avec le fichier BioLip2 au format tabulé."""
    
    df_ecoli = pd.read_csv(ecoli_csv)
    
    COLUMNS = [
        'PDB_ID', 'Receptor_chain', 'Resolution', 'Binding_site_code',
        'Ligand_ID', 'Ligand_chain', 'Ligand_serial',
        'Binding_site_residues_PDB', 'Binding_site_residues_renum',
        'Catalytic_site_PDB', 'Catalytic_site_renum',
        'EC_number', 'GO_terms',
        'Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB',
        'UniProt_ID', 'PubMed_ID', 'Ligand_residue_seq', 'Receptor_sequence'
    ]
    
    df_biolip = pd.read_csv(new_biolip_file, sep='\t', header=None, names=COLUMNS, low_memory=False)
    biolip_pdbs = set(df_biolip['PDB_ID'].astype(str).str.upper().str.strip())
    
    def find_matching_pdbs(pdb_cell):
        if pd.isna(pdb_cell) or pdb_cell == '':
            return []
        pdbs_ecoli = [p.strip().upper() for p in str(pdb_cell).split(';') if len(p.strip()) >= 3]
        return [p for p in pdbs_ecoli if p in biolip_pdbs]
    
    df_ecoli['matching_pdbs'] = df_ecoli['pdb_ids'].apply(find_matching_pdbs)
    df_ecoli['n_matches'] = df_ecoli['matching_pdbs'].apply(len)
    df_ecoli['has_biolip2_match'] = df_ecoli['n_matches'] > 0
    
    df_matches = df_ecoli[df_ecoli['has_biolip2_match']].copy()
    df_matches['pdb_ids_all'] = df_matches['pdb_ids']
    df_matches['pdb_ids'] = df_matches['matching_pdbs'].apply(lambda x: ';'.join(x))
    df_matches_clean = df_matches.drop(columns=['matching_pdbs', 'n_matches', 'has_biolip2_match'], errors='ignore')
    
    df_matches_clean.to_csv(output_csv, index=False)
    return df_matches_clean


df_matched = match_ecoli_with_new_biolip2(
    ecoli_csv="Résultats/Résultats finaux propres/Escherichia_Coli_K12_protéome_complet.csv",
    new_biolip_file="BioLiP.txt",
    output_csv="ecoli_in_new_biolip2.csv"
)


#%% Ecoli biolip2 2026 - ecoli biolip2 article 


def subtract_old_intersection(new_csv, old_intersection_csv, output_csv="ecoli_new_biolip2_only.csv"):
    """Retire du nouveau dataset les protéines déjà présentes dans l'ancienne intersection."""
    
    df_new = pd.read_csv(new_csv)
    df_old = pd.read_csv(old_intersection_csv)
    old_uniprots = set(df_old['uniprot_id'].astype(str).str.upper().str.strip())
    
    df_new['uniprot_upper'] = df_new['uniprot_id'].astype(str).str.upper().str.strip()
    mask_new_only = ~df_new['uniprot_upper'].isin(old_uniprots)
    df_unique = df_new[mask_new_only].copy()
    df_unique = df_unique.drop(columns=['uniprot_upper'], errors='ignore')
    
    df_unique.to_csv(output_csv, index=False)
    return df_unique


df_nouvelles = subtract_old_intersection(
    new_csv="Escherichia_Coli_Protéines_présentes_dans_Biolip2_2026.csv",
    old_intersection_csv="Escherichia_Coli_Protéines_présentes_dans_Biolip2_article.csv",
    output_csv="ecoli_new_biolip2_ONLY.csv"
)


#%% Conversion du fichier Ecoli biolip2 2026 vs article pour avoir les informations de Biolip2


def filter_biolip2_for_ecoli_csv(ecoli_csv, biolip_file, output_file="BioLiP_ecoli_new_only.csv"):
    """Filtre BioLiP2.txt pour ne garder que les entrées E. coli."""
    
    df_ecoli = pd.read_csv(ecoli_csv)
    ecoli_uniprots = set(df_ecoli['uniprot_id'].astype(str).str.upper().str.strip())
    
    all_pdbs = set()
    for pdb_cell in df_ecoli['pdb_ids'].dropna():
        if isinstance(pdb_cell, str):
            pdbs = [p.strip().upper() for p in pdb_cell.split(';') if p.strip()]
            all_pdbs.update(pdbs)
    
    COLUMNS = [
        'PDB_ID', 'Receptor_chain', 'Resolution', 'Binding_site_code',
        'Ligand_ID', 'Ligand_chain', 'Ligand_serial',
        'Binding_site_residues_PDB', 'Binding_site_residues_renum',
        'Catalytic_site_PDB', 'Catalytic_site_renum',
        'EC_number', 'GO_terms',
        'Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB',
        'UniProt_ID', 'PubMed_ID', 'Ligand_residue_seq', 'Receptor_sequence'
    ]
    
    df_biolip = pd.read_csv(biolip_file, sep='\t', header=None, names=COLUMNS, low_memory=False)
    
    df_biolip['PDB_ID_clean'] = df_biolip['PDB_ID'].astype(str).str.upper().str.strip()
    df_biolip['UniProt_ID_clean'] = df_biolip['UniProt_ID'].astype(str).str.upper().str.strip()
    
    mask = (
        df_biolip['PDB_ID_clean'].isin(all_pdbs) |
        df_biolip['UniProt_ID_clean'].isin(ecoli_uniprots)
    )
    df_filtered = df_biolip[mask].copy()
    df_filtered = df_filtered.drop(columns=['PDB_ID_clean', 'UniProt_ID_clean'])
    
    if len(df_filtered) > 0:
        df_filtered.to_csv(output_file, sep=',', index=False, header=True, na_rep='')
    
    return df_filtered


df_biolip_ecoli = filter_biolip2_for_ecoli_csv(
    ecoli_csv="Résultats/Résultats finaux propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Ecoli.csv",
    biolip_file="BioLiP.txt",
    output_file="BioLiP_ecoli_new_only.csv"
)


#%% Compter le nombre de PDB dans le génome


df = pd.read_csv('BioLiP_ecoli_new_only.csv')

total_pdb = df['PDB_ID'].dropna().apply(
    lambda x: len(str(x).split(';')) if str(x).strip() else 0
).sum()

print(f"Nombre total de PDB : {total_pdb}")


#%% Filtrer avec les IC50


import re

INPUT_CSV = "BioLiP_ecoli_new_only.csv"
OUTPUT_CSV = "BioLiP_ecoli_ki_kd_ic50.csv"
AFF_COLS = ['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']


def detect_affinity_type(row):
    for col in AFF_COLS:
        if col not in row or pd.isna(row[col]):
            continue
        val = str(row[col]).upper()
        has_number = bool(re.search(r'\d+\.?\d*\s*(fM|pM|nM|uM|µM|mM)', str(row[col]), re.IGNORECASE))
        if not has_number:
            continue
        if re.search(r'\bKI\b', val) and not re.search(r'\bIC50\b|\bEC50\b', val):
            return 'Ki'
        elif re.search(r'\bKD\b', val) and not re.search(r'\bIC50\b|\bEC50\b', val):
            return 'Kd'
        elif re.search(r'\bIC50\b', val):
            return 'IC50'
    return None


def extract_value_molar(aff_str, target_type):
    if pd.isna(aff_str) or str(aff_str).strip() == '':
        return np.nan
    val_upper = str(aff_str).upper()
    if target_type in ('Ki', 'Kd'):
        pattern = r'\b' + target_type.upper() + r'\b'
        if not re.search(pattern, val_upper):
            return np.nan
        if re.search(r'\bIC50\b|\bEC50\b', val_upper):
            return np.nan
    elif target_type == 'IC50':
        if not re.search(r'\bIC50\b', val_upper):
            return np.nan
    match = re.search(r'(\d+\.?\d*[eE]?[+-]?\d*)\s*(fM|pM|nM|uM|µM|mM)', str(aff_str), re.IGNORECASE)
    if not match:
        return np.nan
    val = float(match.group(1))
    unit = match.group(2).lower()
    conversions = {'fm': 1e-15, 'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'µm': 1e-6, 'mm': 1e-3}
    molar = val * conversions.get(unit, np.nan)
    if pd.isna(molar) or molar <= 0:
        return np.nan
    return molar


df = pd.read_csv(INPUT_CSV, low_memory=False)
df['affinity_type'] = df.apply(detect_affinity_type, axis=1)
df_filtered = df[df['affinity_type'].notna()].copy()

valid_cols = [c for c in AFF_COLS if c in df_filtered.columns]

ki_kd_values = []
for _, row in df_filtered.iterrows():
    best = np.nan
    for col in valid_cols:
        for t in ('Ki', 'Kd'):
            m = extract_value_molar(row[col], t)
            if not pd.isna(m):
                p = -np.log10(m)
                if pd.isna(best) or p > best:
                    best = p
    ki_kd_values.append(best)
df_filtered['pKi_Kd'] = ki_kd_values

ic50_values = []
for _, row in df_filtered.iterrows():
    best = np.nan
    for col in valid_cols:
        m = extract_value_molar(row[col], 'IC50')
        if not pd.isna(m):
            p = -np.log10(m)
            if pd.isna(best) or p > best:
                best = p
    ic50_values.append(best)
df_filtered['pIC50'] = ic50_values

df_filtered.to_csv(OUTPUT_CSV, index=False)


#%% observation des données


import matplotlib.pyplot as plt

df = pd.read_csv(INPUT_CSV, low_memory=False)
df['affinity_type'] = df.apply(detect_affinity_type, axis=1)
df_filtered = df[df['affinity_type'].notna()].copy()

valid_cols = [c for c in AFF_COLS if c in df_filtered.columns]

ki_kd_values = []
for _, row in df_filtered.iterrows():
    best = np.nan
    for col in valid_cols:
        for t in ('Ki', 'Kd'):
            m = extract_value_molar(row[col], t)
            if not pd.isna(m):
                p = -np.log10(m)
                if pd.isna(best) or p > best:
                    best = p
    ki_kd_values.append(best)
df_filtered['pKi_Kd'] = ki_kd_values

ic50_values = []
for _, row in df_filtered.iterrows():
    best = np.nan
    for col in valid_cols:
        m = extract_value_molar(row[col], 'IC50')
        if not pd.isna(m):
            p = -np.log10(m)
            if pd.isna(best) or p > best:
                best = p
    ic50_values.append(best)
df_filtered['pIC50'] = ic50_values

df_filtered.to_csv(OUTPUT_CSV, index=False)

pki_kd = df_filtered['pKi_Kd'].dropna()
pic50 = df_filtered['pIC50'].dropna()


def classify_potency(series, name):
    bins = {
        'Très puissant (< 1 nM, pAff ≥ 9)': (series >= 9).sum(),
        'Puissant (1-100 nM, pAff 7-9)': ((series >= 7) & (series < 9)).sum(),
        'Modéré (0.1-1 µM, pAff 6-7)': ((series >= 6) & (series < 7)).sum(),
        'Faible (1-10 µM, pAff 5-6)': ((series >= 5) & (series < 6)).sum(),
        'Très faible (> 10 µM, pAff < 5)': (series < 5).sum(),
    }
    total = len(series)
    for label, count in bins.items():
        pct = count / total * 100
        bar = '█' * int(pct / 2)
        print(f"   {label:45s} : {count:4d} ({pct:5.1f}%) {bar}")
    print(f"\n   Médiane : {series.median():.2f}")
    print(f"   Moyenne : {series.mean():.2f}")
    print(f"   Std     : {series.std():.2f}")
    return bins


bins_ki_kd = classify_potency(pki_kd, "pKi/Kd") if len(pki_kd) > 0 else None
bins_ic50 = classify_potency(pic50, "pIC50") if len(pic50) > 0 else None

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Distribution des affinités — BioLiP2 E. coli", fontsize=14, fontweight='bold')

color_ki = '#2196F3'
color_ic50 = '#FF9800'

ax1 = axes[0, 0]
if len(pki_kd) > 0:
    ax1.hist(pki_kd, bins=30, density=True, alpha=0.7, color=color_ki, edgecolor='white')
    ax1.axvline(pki_kd.median(), color='red', linestyle='--', linewidth=2)
    ax1.set_xlabel('pKi/Kd')
    ax1.set_ylabel('Densité')
    ax1.set_title(f'pKi/Kd (n={len(pki_kd)})')

ax2 = axes[0, 1]
if len(pic50) > 0:
    ax2.hist(pic50, bins=30, density=True, alpha=0.7, color=color_ic50, edgecolor='white')
    ax2.axvline(pic50.median(), color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('pIC50')
    ax2.set_ylabel('Densité')
    ax2.set_title(f'pIC50 (n={len(pic50)})')

ax3 = axes[1, 0]
if len(pki_kd) > 0:
    ax3.hist(pki_kd, bins=30, density=True, alpha=0.5, color=color_ki, edgecolor='white', label=f'pKi/Kd')
if len(pic50) > 0:
    ax3.hist(pic50, bins=30, density=True, alpha=0.5, color=color_ic50, edgecolor='white', label=f'pIC50')
ax3.set_xlabel('pAffinity')
ax3.set_ylabel('Densité')
ax3.set_title('Comparaison pKi/Kd vs pIC50')
ax3.legend()

ax4 = axes[1, 1]
categories = ['< 1 nM', '1-100 nM', '0.1-1 µM', '1-10 µM', '> 100 µM']
x = np.arange(len(categories))
width = 0.35

if bins_ki_kd and len(pki_kd) > 0:
    vals_ki = [v / len(pki_kd) * 100 for v in bins_ki_kd.values()]
    ax4.bar(x - width/2, vals_ki, width, label='pKi/Kd', color=color_ki, alpha=0.7)
if bins_ic50 and len(pic50) > 0:
    vals_ic50 = [v / len(pic50) * 100 for v in bins_ic50.values()]
    ax4.bar(x + width/2, vals_ic50, width, label='pIC50', color=color_ic50, alpha=0.7)

ax4.set_xlabel('Puissance')
ax4.set_ylabel('Proportion (%)')
ax4.set_title('Répartition par catégorie')
ax4.set_xticks(x)
ax4.set_xticklabels(categories, fontsize=8)
ax4.legend()

plt.tight_layout()
plt.savefig('distribution_affinites_ecoli.png', dpi=300, bbox_inches='tight')
plt.show()


#%% Vérification cofacteurs


df = pd.read_csv("tier2_ecoli_inhibition_set.csv")

COFACTORS = {
    'NAD', 'NAP', 'FAD', 'FMN', 'COA', 'SAM', 'SAH',
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',
    'HEM', 'HEC', 'PLP', 'TPP', 'B12',
    'MG', 'ZN', 'CA', 'FE', 'MN', 'CU',
    'SO4', 'PO4', 'GOL', 'EDO', 'DMS', 'ACT'
}

if 'Ligand_ID' in df.columns:
    df['is_cofactor'] = df['Ligand_ID'].str.upper().isin(COFACTORS)
    df_clean = df[~df['is_cofactor']].copy()
    print(f"Paires restantes après retrait cofacteurs : {len(df_clean)}")


#%% Séparation du dataset E coli Ki Kd IC50 en 3 datasets


INPUT_CSV = "BioLiP_ecoli_ki_kd_ic50.csv"
df = pd.read_csv(INPUT_CSV, low_memory=False)

tier2 = df[df['affinity_type'] == 'IC50'].copy()


def classify_inhibitor(pic50):
    if pd.isna(pic50):
        return 'Non parsé'
    elif pic50 >= 8:
        return 'Hit avancé / Lead (< 10 nM)'
    elif pic50 >= 7:
        return 'Hit puissant (10-100 nM)'
    elif pic50 >= 6:
        return 'Hit confirmé (0.1-1 µM)'
    elif pic50 >= 5:
        return 'Hit précoce (1-10 µM)'
    elif pic50 >= 4:
        return 'Hit faible (10-100 µM)'
    else:
        return 'Inactif (> 100 µM)'


tier2['inhibitor_class'] = tier2['pIC50'].apply(classify_inhibitor)

target_summary = tier2.groupby(['UniProt_ID', 'PDB_ID']).agg(
    n_inhibitors=('Ligand_ID', 'nunique'),
    best_pIC50=('pIC50', 'max'),
    median_pIC50=('pIC50', 'median')
).reset_index().sort_values('best_pIC50', ascending=False)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('Tier 2 : E. coli K12 Inhibition Set — pIC50', fontsize=13, fontweight='bold')

ax1 = axes[0]
pic50_valid = tier2['pIC50'].dropna()
if len(pic50_valid) > 0:
    ax1.hist(pic50_valid, bins=20, density=True, alpha=0.7, color='#FF5722', edgecolor='white')
    ax1.axvline(pic50_valid.median(), color='black', linestyle='--', linewidth=2)
    ax1.set_xlabel('pIC50')
    ax1.set_ylabel('Densité')
    ax1.set_title(f'Distribution (n={len(pic50_valid)})')

ax2 = axes[1]
class_order = ['Hit avancé / Lead (< 10 nM)', 'Hit puissant (10-100 nM)',
               'Hit confirmé (0.1-1 µM)', 'Hit précoce (1-10 µM)',
               'Hit faible (10-100 µM)', 'Inactif (> 100 µM)']
short_labels = ['< 10 nM', '10-100 nM', '0.1-1 µM', '1-10 µM', '10-100 µM', '> 100 µM']
colors_bar = ['#1B5E20', '#4CAF50', '#8BC34A', '#FFC107', '#FF9800', '#F44336']

counts = [tier2[tier2['inhibitor_class'] == c].shape[0] for c in class_order]
bars = ax2.bar(range(len(class_order)), counts, color=colors_bar, edgecolor='white')
ax2.set_xticks(range(len(class_order)))
ax2.set_xticklabels(short_labels, fontsize=7)
ax2.set_ylabel('Nombre de paires')
ax2.set_title('Classification pharmacologique')
for bar, count in zip(bars, counts):
    if count > 0:
        ax2.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + 0.5,
                 str(count), ha='center', va='bottom', fontsize=9)

ax3 = axes[2]
top_targets = target_summary.head(10)
if len(top_targets) > 0:
    y_pos = range(len(top_targets))
    ax3.barh(y_pos, top_targets['best_pIC50'].values, color='#FF5722', alpha=0.7)
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels([f"{r['UniProt_ID']}\n({r['PDB_ID']})" for _, r in top_targets.iterrows()], fontsize=7)
    ax3.set_xlabel('Meilleur pIC50')
    ax3.set_title('Top 10 cibles')

plt.tight_layout()
plt.savefig('tier2_ecoli_inhibition_set.png', dpi=300, bbox_inches='tight')
plt.show()

tier1 = df[df['affinity_type'].isin(['Ki', 'Kd'])].copy()
tier3 = df.copy()

tier1.to_csv('tier1_ecoli_gold_standard.csv', index=False)
tier2.to_csv('tier2_ecoli_inhibition_set.csv', index=False)
tier3.to_csv('tier3_ecoli_full.csv', index=False)

print(f"Tier 1 : {len(tier1)} lignes")
print(f"Tier 2 : {len(tier2)} lignes")
print(f"Tier 3 : {len(tier3)} lignes")


#%% Enlever les doublons de PDB


df = pd.read_csv("Résultats/Résultats finaux propres/BioLiP_ecoli_ki_kd_ic50.csv")

print(f"Avant déduplication : {len(df)} lignes")

df["Resolution"] = pd.to_numeric(df["Resolution"], errors="coerce")
df_dedup = (df
    .sort_values("Resolution", ascending=True)
    .drop_duplicates(subset=["UniProt_ID", "Ligand_ID"], keep="first")
)

print(f"Après déduplication : {len(df_dedup)} lignes")
print(f"Retirés : {len(df) - len(df_dedup)}")

df_dedup.to_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv", index=False)


#%% Compter les Ki kd IC50


df = pd.read_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv")

n_ki_kd = df["pKi_Kd"].notna().sum()
n_ic50 = df["pIC50"].notna().sum()
n_both = (df["pKi_Kd"].notna() & df["pIC50"].notna()).sum()
n_none = (df["pKi_Kd"].isna() & df["pIC50"].isna()).sum()

print(f"Total      : {len(df)}")
print(f"Ki/Kd      : {n_ki_kd}")
print(f"IC50       : {n_ic50}")
print(f"Les deux   : {n_both}")
print(f"Sans rien  : {n_none}")


#%% Nettoie doublons dans fichier basé sur PDB


import os

INPUT_CSV = "Code/Récupération_base_de_données_Chembl/Résultats/Résultats_finaux_propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Biolip2.csv"
OUTPUT_CSV = "Code/Récupération_base_de_données_Chembl/Résultats/Résultats_finaux_propres/Escherichia_coli_biolip2_2026_vs_biolip2_article_format_Biolip2_CLEANED.csv"
PDB_COLUMN = "PDB_ID"
PROTEIN_COLUMN = "UniProt_ID"


def clean_duplicates():
    df = pd.read_csv(INPUT_CSV, dtype=str, low_memory=False)
    
    if PDB_COLUMN not in df.columns:
        print(f"Colonne '{PDB_COLUMN}' introuvable.")
        print(f"Colonnes : {df.columns.tolist()}")
        return
    
    df[PDB_COLUMN] = df[PDB_COLUMN].str.strip().str.upper()
    df = df[df[PDB_COLUMN] != "NAN"]
    df = df[df[PDB_COLUMN] != ""]
    
    initial_count = len(df)
    df_cleaned = df.drop_duplicates(subset=[PDB_COLUMN], keep='first')
    removed_count = initial_count - len(df_cleaned)
    
    print(f"Doublons supprimés : {removed_count}")
    print(f"Lignes finales : {len(df_cleaned)}")
    
    output_dir = os.path.dirname(OUTPUT_CSV)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    df_cleaned.to_csv(OUTPUT_CSV, index=False)
    print(f"Sauvegardé : {OUTPUT_CSV}")


if __name__ == "__main__":
    clean_duplicates()