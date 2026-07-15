import pandas as pd

NEW_BIOLIP_FILE = "BioLiP.txt"
OLD_CSV_FILE = "pdb_extracted.csv"
OUTPUT_FILE = "Biolip_new_vs_Biolp_article.csv"

COLUMNS = [
    'PDB_ID', 'Receptor_chain', 'Resolution', 'Binding_site_code',
    'Ligand_ID', 'Ligand_chain', 'Ligand_serial',
    'Binding_site_residues_PDB', 'Binding_site_residues_renum',
    'Catalytic_site_PDB', 'Catalytic_site_renum',
    'EC_number', 'GO_terms',
    'Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB',
    'UniProt_ID', 'PubMed_ID', 'Ligand_residue_seq', 'Receptor_sequence'
]

df_new = pd.read_csv(NEW_BIOLIP_FILE, sep='\t', header=None, names=COLUMNS, low_memory=False)
df_old = pd.read_csv(OLD_CSV_FILE)

pdb_col_old = None
for col in df_old.columns:
    if 'pdb' in col.lower():
        pdb_col_old = col
        break

if not pdb_col_old:
    raise ValueError(f"Colonne PDB non trouvée. Colonnes disponibles: {list(df_old.columns)}")

df_new['PDB_ID_CLEAN'] = df_new['PDB_ID'].str.upper().str.strip()
df_old['PDB_ID_CLEAN'] = df_old[pdb_col_old].astype(str).str.upper().str.strip()

old_pdbs = set(df_old['PDB_ID_CLEAN'].unique())
new_pdbs = set(df_new['PDB_ID_CLEAN'].unique())
nouveaux_pdbs = new_pdbs - old_pdbs

print(f"Nouveaux PDBs : {len(nouveaux_pdbs)}")

if len(nouveaux_pdbs) > 0:
    df_nouveaux = df_new[df_new['PDB_ID_CLEAN'].isin(nouveaux_pdbs)].copy()
    df_nouveaux = df_nouveaux.drop(columns=['PDB_ID_CLEAN'])
    print(f"{len(df_nouveaux)} interactions pour {len(nouveaux_pdbs)} nouveaux PDBs")
    print(f"Exemples : {list(nouveaux_pdbs)[:5]}")


#%% Application de filtres comme l'article

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
        match = re.search(r'(\d+\.?\d*)\s*(nM|uM|µM|mM|pM|fM)', str(aff_str), re.IGNORECASE)
        if not match:
            return np.nan
        val, unit = float(match.group(1)), match.group(2).lower()
        conversions = {'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'µm': 1e-6, 'mm': 1e-3}
        molar = val * conversions.get(unit, 1e-9)
        return -np.log10(molar) if molar > 0 else np.nan

    aff_cols = ['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']
    for col in aff_cols:
        df[f'{col}_pKi'] = df[col].apply(extract_paffinity)

    df['pAffinity'] = df[[f'{c}_pKi' for c in aff_cols]].max(axis=1, skipna=True)
    df['pAffinity_rank'] = df['pAffinity'].fillna(-999)

    print(f"Affinités parsées : {df['pAffinity'].notna().sum()}/{len(df)}")
    return df


def cluster_sequences_cdhit(df, sequence_col='Receptor_sequence', threshold=0.9):
    """Clusterise les séquences à threshold% d'identité avec CD-HIT."""
    try:
        subprocess.run(['cd-hit', '-h'], capture_output=True, text=True)
        cdhit_available = True
    except FileNotFoundError:
        cdhit_available = False
        print("CD-HIT non trouvé. Fallback sur UniPROT ID.")

    records = []
    valid_indices = []

    for idx, row in df.iterrows():
        seq = str(row[sequence_col]).strip()
        if len(seq) > 20:
            rec = SeqRecord(
                Seq(seq),
                id=f"idx{idx}|{row.get('UniProt_ID', 'unknown')}",
                description=""
            )
            records.append(rec)
            valid_indices.append(idx)

    if len(records) == 0:
        return {i: [i] for i in df.index}

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_in = os.path.join(tmpdir, "input.fasta")
        fasta_out = os.path.join(tmpdir, "output")
        SeqIO.write(records, fasta_in, "fasta")

        if cdhit_available:
            cmd = [
                'cd-hit', '-i', fasta_in, '-o', fasta_out,
                '-c', str(threshold), '-n', '5', '-d', '0', '-s', '0.9', '-S', '0'
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                return {i: [i] for i in valid_indices}
            return parse_cdhit_clusters(fasta_out + ".clstr", df.index)
        else:
            groups = df.groupby('UniProt_ID').apply(lambda x: list(x.index)).to_dict()
            return {f"cluster_{k}": v for k, v in groups.items()}


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
                match = re.search(r'>idx(\d+)\|', line)
                if match and current_cluster:
                    idx = int(match.group(1))
                    if current_cluster:
                        clusters[current_cluster].append(idx)

    valid_clusters = {k: [i for i in v if i in all_indices]
                      for k, v in clusters.items() if v}
    print(f"Clustering terminé : {len(valid_clusters)} clusters")
    return valid_clusters


def filter_protein_redundancy(df, clusters):
    """Pour chaque cluster, garde l'entrée avec la meilleure affinité."""
    indices_to_keep = []

    for cluster_id, idx_list in clusters.items():
        if len(idx_list) == 1:
            indices_to_keep.extend(idx_list)
        else:
            sub_df = df.loc[idx_list]
            best_idx = sub_df['pAffinity_rank'].idxmax()
            indices_to_keep.append(best_idx)

    result = df.loc[indices_to_keep].copy()
    print(f"Filtre redondance : {len(df)} -> {len(result)} ({len(df) - len(result)} supprimées)")
    return result


def main_pipeline(input_csv, output_csv):
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"Entrées initiales : {len(df)}")

    df = parse_affinities(df)
    clusters = cluster_sequences_cdhit(df, threshold=0.9)
    df_final = filter_protein_redundancy(df, clusters)

    cols_to_drop = [c for c in df_final.columns if c.endswith('_pKi') or c in ['pAffinity_rank']]
    df_final = df_final.drop(columns=cols_to_drop, errors='ignore')
    df_final.to_csv(output_csv, index=False)

    strong = (df_final['pAffinity'] >= 6).sum()
    print(f"Résultat : {len(df_final)} paires, dont {strong} avec pKi >= 6")


if __name__ == "__main__":
    INPUT = "Résultats/Résultats finaux propres/Biolip_new_vs_Biolp_article.csv"
    OUTPUT = "biolip_filtered_90identity.csv"
    main_pipeline(INPUT, OUTPUT)


#%% Observation des affinités

df = pd.read_csv("Résultats/Résultats finaux propres/Biolip2_2026_vs_biolip2_article_plusieurs_ligands.csv")

total = len(df)
with_affinity = df['pAffinity'].notna().sum()
without_affinity = df['pAffinity'].isna().sum()

print(f"Total paires     : {total}")
print(f"Avec affinité    : {with_affinity} ({with_affinity / total * 100:.1f}%)")
print(f"Sans affinité    : {without_affinity} ({without_affinity / total * 100:.1f}%)")

if with_affinity > 0:
    print(f"pKi >= 8  : {(df['pAffinity'] >= 8).sum()}")
    print(f"6 <= pKi < 8 : {((df['pAffinity'] >= 6) & (df['pAffinity'] < 8)).sum()}")
    print(f"5 <= pKi < 6 : {((df['pAffinity'] >= 5) & (df['pAffinity'] < 6)).sum()}")
    print(f"pKi < 5   : {(df['pAffinity'] < 5).sum()}")


#%% Vérification cofacteurs

df = pd.read_csv("biolip_filtered_90identity.csv")

COFACTORS = {
    'NAD', 'NAP', 'FAD', 'FMN', 'COA', 'SAM', 'SAH',
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',
    'HEM', 'HEC', 'PLP', 'TPP', 'B12',
    'MG', 'ZN', 'CA', 'FE', 'MN', 'CU',
    'SO4', 'PO4', 'GOL', 'EDO', 'DMS', 'ACT'
}

if 'Ligand_ID' in df.columns:
    df['is_cofactor'] = df['Ligand_ID'].str.upper().isin(COFACTORS)
    n_cofactors = df['is_cofactor'].sum()
    print(f"Cofacteurs détectés : {n_cofactors}/{len(df)}")

    if n_cofactors > 0:
        print(df[df['is_cofactor']]['Ligand_ID'].value_counts())

    df_clean = df[~df['is_cofactor']].copy()
    print(f"Après nettoyage : {len(df_clean)} paires")
    df_clean.to_csv("Biolip2_2026_vs_biolip2_article_plusieurs_ligands.csv", index=False)


#%% Version des filtres comme l'article avec filtre Ki/Kd avant


def filter_ki_kd_strict(df, aff_cols=['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']):
    """Garde uniquement les lignes avec Ki ou Kd explicites (pas IC50/EC50)."""

    def has_ki_or_kd(row):
        for col in aff_cols:
            if col not in row or pd.isna(row[col]):
                continue
            val = str(row[col]).upper()
            has_ki_kd = bool(re.search(r'\bKI\b|\bKD\b', val))
            has_ic50_ec50 = bool(re.search(r'\bIC50\b|\bEC50\b', val))
            if has_ki_kd and not has_ic50_ec50:
                if re.search(r'\d+\.?\d*\s*(pM|nM|µM|uM|mM|fM)', val, re.IGNORECASE):
                    return True
        return False

    mask = df.apply(has_ki_or_kd, axis=1)
    df_filtered = df[mask].copy()
    print(f"Filtre Ki/Kd : {len(df)} -> {len(df_filtered)} ({len(df) - len(df_filtered)} retirés)")
    return df_filtered


def parse_affinities(df, aff_cols=['Affinity_manual', 'Affinity_MOAD', 'Affinity_PDBbind', 'Affinity_BindingDB']):
    """Convertit les affinités Ki/Kd en pKi/pKd."""

    def extract_paffinity(aff_str):
        if pd.isna(aff_str) or aff_str == '' or str(aff_str).upper() == 'NAN':
            return np.nan
        val_upper = str(aff_str).upper()
        if not (re.search(r'\bKI\b|\bKD\b', val_upper) and not re.search(r'\bIC50\b|\bEC50\b', val_upper)):
            return np.nan
        match = re.search(r'(\d+\.?\d*[eE]?[+-]?\d*)\s*(pM|nM|uM|µM|mM|fM)', str(aff_str), re.IGNORECASE)
        if not match:
            return np.nan
        val, unit = float(match.group(1)), match.group(2).lower()
        conversions = {'pm': 1e-12, 'nm': 1e-9, 'um': 1e-6, 'µm': 1e-6, 'mm': 1e-3, 'fm': 1e-15}
        molar = val * conversions.get(unit, np.nan)
        if pd.isna(molar) or molar <= 0:
            return np.nan
        return -np.log10(molar)

    valid_cols = [c for c in aff_cols if c in df.columns]
    for col in valid_cols:
        df[f'{col}_pKi'] = df[col].apply(extract_paffinity)

    df['pAffinity'] = df[[f'{c}_pKi' for c in valid_cols]].max(axis=1, skipna=True)
    df['pAffinity_rank'] = df['pAffinity'].fillna(-999)

    valid = df['pAffinity'].notna().sum()
    print(f"Affinités parsées : {valid}/{len(df)}")
    return df


def cluster_sequences_cdhit(df, sequence_col='Receptor_sequence', threshold=0.9):
    """Clusterise les séquences à threshold% d'identité avec CD-HIT."""
    try:
        subprocess.run(['cd-hit', '-h'], capture_output=True, text=True)
        cdhit_available = True
    except FileNotFoundError:
        cdhit_available = False
        print("CD-HIT non trouvé. Fallback sur séquence exacte.")

    if sequence_col not in df.columns:
        seq_candidates = [c for c in df.columns if 'seq' in c.lower()]
        if seq_candidates:
            sequence_col = seq_candidates[0]

    records = []
    valid_indices = []

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
    print(f"Clustering terminé : {len(valid_clusters)} clusters")
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
    print(f"Filtre redondance : {len(df)} -> {len(result)} ({len(df) - len(result)} doublons)")
    return result


def main_pipeline(input_csv, output_csv):
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"Entrées initiales : {len(df)}")

    df = filter_ki_kd_strict(df)
    if len(df) == 0:
        print("Aucune donnée Ki/Kd exploitable.")
        return

    df = parse_affinities(df)
    if df['pAffinity'].notna().sum() == 0:
        print("Aucune affinité Ki/Kd parsable.")
        return

    clusters = cluster_sequences_cdhit(df, threshold=0.9)
    df_final = filter_protein_redundancy(df, clusters)

    cols_to_drop = [c for c in df_final.columns if c.endswith('_pKi') or c == 'pAffinity_rank']
    df_export = df_final.drop(columns=cols_to_drop, errors='ignore')
    df_export.to_csv(output_csv, index=False)

    print(f"Résultat final : {len(df_final)} paires, {(df_final['pAffinity'] >= 6).sum()} avec pKi >= 6")
    print(f"Réduction : {len(df)} -> {len(df_final)} ({(1 - len(df_final) / len(df)) * 100:.1f}%)")


if __name__ == "__main__":
    INPUT = "Résultats/Résultats finaux propres/Biolip_new_vs_Biolp_article.csv"
    OUTPUT = "biolip_filtered_90identity.csv"
    main_pipeline(INPUT, OUTPUT)


#%% Filtre anti leakage (30% Idt avec biolip2 article)

import urllib.request
import gzip
import shutil
import time
import requests

VALIDATION_FILE = "Résultats/Résultats finaux propres/BioLiP_ecoli_ki_kd_ic50.csv"
FINETUNE_FILE = "PDB finetuning article.csv"
FINETUNE_PDB_COL = "PDB_ID"
OUTPUT_DIR = "leakage_check_E.coli"
IDENTITY_THRESHOLD = 0.3

SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
SIFTS_GZ = "pdb_chain_uniprot.tsv.gz"
SIFTS_TSV = "pdb_chain_uniprot.tsv"


def get_sifts():
    if not os.path.exists(SIFTS_TSV):
        if not os.path.exists(SIFTS_GZ):
            print("Téléchargement SIFTS (~80 Mo)...")
            urllib.request.urlretrieve(SIFTS_URL, SIFTS_GZ)
        with gzip.open(SIFTS_GZ, "rb") as f_in, open(SIFTS_TSV, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

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
    return df[["PDB", "SP_PRIMARY"]].drop_duplicates()


def load_validation(filepath):
    df = pd.read_csv(filepath)
    df["_pdb"] = df["PDB_ID"].astype(str).str.strip().str.lower()
    df["_uniprot"] = df["UniProt_ID"].astype(str).str.strip()
    pdbs = set(df["_pdb"].dropna().unique()) - {"nan", ""}
    uniprots = set(df["_uniprot"].dropna().unique()) - {"nan", ""}
    print(f"Validation : {len(df)} lignes, {len(pdbs)} PDB, {len(uniprots)} UniProt")
    return df, pdbs, uniprots


def load_finetune(filepath, sifts):
    df = pd.read_csv(filepath)
    df["_pdb"] = df[FINETUNE_PDB_COL].astype(str).str.strip().str.lower()
    pdbs = set(df["_pdb"].dropna().unique()) - {"nan", ""}
    merged = df.merge(sifts, left_on="_pdb", right_on="PDB", how="left")
    uniprots = set(merged["SP_PRIMARY"].dropna().unique())
    print(f"Finetuning : {len(pdbs)} PDB -> {len(uniprots)} UniProt")
    return df, pdbs, uniprots


def check_exact_overlap(val_pdbs, val_uniprots, ft_pdbs, ft_uniprots):
    pdb_overlap = val_pdbs & ft_pdbs
    uni_overlap = val_uniprots & ft_uniprots

    print(f"PDB en commun     : {len(pdb_overlap)}")
    print(f"UniProt en commun : {len(uni_overlap)}")

    if pdb_overlap:
        for p in sorted(list(pdb_overlap)[:20]):
            print(f"  {p}")
    if uni_overlap:
        for u in sorted(list(uni_overlap)[:30]):
            print(f"  {u}")

    return pdb_overlap, uni_overlap


def download_sequences(uniprot_ids, output_fasta, label):
    print(f"Téléchargement séquences {label} ({len(uniprot_ids)})...")

    if os.path.exists(output_fasta):
        n = sum(1 for l in open(output_fasta) if l.startswith(">"))
        if n >= len(uniprot_ids) * 0.8:
            print(f"  {n} séquences déjà présentes, skip.")
            return
        os.remove(output_fasta)

    ids = sorted(uniprot_ids)
    n_ok = 0
    n_fail = 0

    with open(output_fasta, "w") as fout:
        for i in range(0, len(ids), 50):
            batch = ids[i:i + 50]
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
                print(f"  Erreur batch {i}: {e}")
                n_fail += len(batch)

            done = min(i + 50, len(ids))
            if done % 500 == 0 or done == len(ids):
                print(f"  {done}/{len(ids)} (ok: {n_ok}, fail: {n_fail})")
            time.sleep(0.2)

    print(f"  -> {n_ok} séquences dans {output_fasta}")


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


def check_sequence_identity(val_seqs, ft_seqs, threshold=0.3):
    K = 3
    kmer_threshold = threshold * 0.7

    ft_kmers = {}
    for acc, seq in ft_seqs.items():
        if len(seq) >= K:
            ft_kmers[acc] = set(seq[i:i + K] for i in range(len(seq) - K + 1))

    matches = {}
    t0 = time.time()

    for idx, (val_acc, val_seq) in enumerate(val_seqs.items()):
        if len(val_seq) < K:
            continue

        val_kset = set(val_seq[i:i + K] for i in range(len(val_seq) - K + 1))
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
            print(f"  {idx + 1}/{len(val_seqs)} leaks={len(matches)} ({elapsed:.0f}s)")

    n_clean = len(val_seqs) - len(matches)
    print(f"Homologues >= 30% : {len(matches)}")
    print(f"Séquences propres : {n_clean}")

    if matches:
        sorted_matches = sorted(matches.items(), key=lambda x: -x[1][1])
        for val_acc, (ft_acc, score) in sorted_matches[:30]:
            print(f"  {val_acc} <-> {ft_acc} (score: {score:.2f})")

    return matches


def export_results(pdb_overlap, uni_overlap, seq_matches, val_df):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    p1 = os.path.join(OUTPUT_DIR, "leakage_summary.txt")
    with open(p1, "w") as f:
        f.write(f"PDB partagés         : {len(pdb_overlap)}\n")
        f.write(f"UniProt partagés     : {len(uni_overlap)}\n")
        f.write(f"Homologues >= 30%    : {len(seq_matches)}\n")
        f.write(f"Total problématiques : {len(uni_overlap | set(seq_matches.keys()))}\n")

    all_leaky = uni_overlap | set(seq_matches.keys())
    p2 = os.path.join(OUTPUT_DIR, "leaky_uniprots.txt")
    with open(p2, "w") as f:
        f.write("uniprot_id\tleak_type\tbest_match\tscore\n")
        for u in sorted(all_leaky):
            if u in uni_overlap and u in seq_matches:
                mt, sc = seq_matches[u]
                f.write(f"{u}\texact+sequence\t{mt}\t{sc:.3f}\n")
            elif u in uni_overlap:
                f.write(f"{u}\texact_uniprot\t-\t1.0\n")
            else:
                mt, sc = seq_matches[u]
                f.write(f"{u}\tsequence_homolog\t{mt}\t{sc:.3f}\n")

    val_clean = val_df[~val_df["_uniprot"].isin(all_leaky)].copy()
    p3 = os.path.join(OUTPUT_DIR, "validation_clean.csv")
    cols_to_drop = [c for c in val_clean.columns if c.startswith("_")]
    val_clean.drop(columns=cols_to_drop, inplace=True)
    val_clean.to_csv(p3, index=False)

    val_leaky = val_df[val_df["_uniprot"].isin(all_leaky)].copy()
    p4 = os.path.join(OUTPUT_DIR, "validation_leaky.csv")
    cols_to_drop = [c for c in val_leaky.columns if c.startswith("_")]
    val_leaky.drop(columns=cols_to_drop, inplace=True)
    val_leaky.to_csv(p4, index=False)

    print(f"{len(all_leaky)} protéines problématiques -> {p2}")
    print(f"{len(val_clean)} entrées propres -> {p3}")
    print(f"{len(val_leaky)} entrées à retirer -> {p4}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    val_fasta = os.path.join(OUTPUT_DIR, "validation_seqs.fasta")
    ft_fasta = os.path.join(OUTPUT_DIR, "finetune_seqs.fasta")

    sifts = get_sifts()
    val_df, val_pdbs, val_uniprots = load_validation(VALIDATION_FILE)
    ft_df, ft_pdbs, ft_uniprots = load_finetune(FINETUNE_FILE, sifts)
    pdb_overlap, uni_overlap = check_exact_overlap(val_pdbs, val_uniprots, ft_pdbs, ft_uniprots)

    val_to_check = val_uniprots - uni_overlap
    download_sequences(val_to_check, val_fasta, "validation")
    download_sequences(ft_uniprots, ft_fasta, "finetuning")

    val_seqs = parse_fasta(val_fasta)
    ft_seqs = parse_fasta(ft_fasta)
    print(f"Séquences parsées : validation={len(val_seqs)}, finetuning={len(ft_seqs)}")

    if len(val_seqs) == 0 or len(ft_seqs) == 0:
        print("Pas assez de séquences téléchargées, check séquence incomplet.")
        seq_matches = {}
    else:
        seq_matches = check_sequence_identity(val_seqs, ft_seqs, IDENTITY_THRESHOLD)

    all_leaky = uni_overlap | set(seq_matches.keys())
    n_clean = len(val_uniprots) - len(all_leaky)

    print(f"Total validation    : {len(val_uniprots)}")
    print(f"PDB identiques      : {len(pdb_overlap)}")
    print(f"UniProt identiques  : {len(uni_overlap)}")
    print(f"Homologues >= 30%   : {len(seq_matches)}")
    print(f"A retirer           : {len(all_leaky)}")
    print(f"Propres             : {n_clean} ({100 * n_clean / max(len(val_uniprots), 1):.1f}%)")

    export_results(pdb_overlap, uni_overlap, seq_matches, val_df)


if __name__ == "__main__":
    main()


#%% Visualisation données ki et IC50

df = pd.read_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv")

PKI_COL = 'pKi_Kd'
IC50_COL = 'Affinity_BindingDB'
UNITE = 'nM'

seuils_pki = [5, 6, 7, 8, 9]

print("--- pKi ---")
for seuil in seuils_pki:
    count = (df[PKI_COL] > seuil).sum()
    print(f"  pKi > {seuil} : {count} ({count / len(df) * 100:.1f}%)")

bins_pki = [0, 5, 6, 7, 8, 9, float('inf')]
labels_pki = ['<=5', '5-6', '6-7', '7-8', '8-9', '>9']
df['pKi_tranche'] = pd.cut(df[PKI_COL], bins=bins_pki, labels=labels_pki, right=True)
print(df['pKi_tranche'].value_counts().sort_index())

print("\n--- IC50 ---")
if UNITE == 'nM':
    for seuil in [1, 10, 100, 1000, 10000]:
        count = (df[IC50_COL] < seuil).sum()
        print(f"  IC50 < {seuil} nM : {count} ({count / len(df) * 100:.1f}%)")
    df['pIC50'] = -np.log10(df[IC50_COL] * 1e-9)
elif UNITE == 'uM':
    for seuil in [0.001, 0.01, 0.1, 1, 10]:
        count = (df[IC50_COL] < seuil).sum()
        print(f"  IC50 < {seuil} µM : {count} ({count / len(df) * 100:.1f}%)")
    df['pIC50'] = -np.log10(df[IC50_COL] * 1e-6)

for seuil in seuils_pki:
    count = (df['pIC50'] > seuil).sum()
    print(f"  pIC50 > {seuil} : {count} ({count / len(df) * 100:.1f}%)")

df.to_csv("resultats_analyse.csv", index=False)


#%% Idem

df = pd.read_csv("BioLiP_ecoli_ki_kd_ic50_sans_doublons.csv")
print(df.columns.tolist())

PIC50_COL = 'pIC50'

df[PIC50_COL] = pd.to_numeric(df[PIC50_COL], errors='coerce')
df_clean = df.dropna(subset=[PIC50_COL]).copy()
print(f"Valeurs valides : {len(df_clean)} / {len(df)}")

seuils = [5, 6, 7, 8, 9]
for seuil in seuils:
    count = (df_clean[PIC50_COL] > seuil).sum()
    print(f"  pIC50 > {seuil} : {count} ({count / len(df_clean) * 100:.1f}%)")

bins = [0, 5, 6, 7, 8, 9, float('inf')]
labels = ['<=5', '5-6', '6-7', '7-8', '8-9', '>9']
df_clean['tranche'] = pd.cut(df_clean[PIC50_COL], bins=bins, labels=labels, right=True)
print(df_clean['tranche'].value_counts().sort_index())

print(f"Min    : {df_clean[PIC50_COL].min():.2f}")
print(f"Max    : {df_clean[PIC50_COL].max():.2f}")
print(f"Moyenne: {df_clean[PIC50_COL].mean():.2f}")
print(f"Médiane: {df_clean[PIC50_COL].median():.2f}")

top20 = df_clean.nlargest(20, PIC50_COL)
print(top20.to_string(index=False))