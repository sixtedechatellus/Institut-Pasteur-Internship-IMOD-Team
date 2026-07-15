# Obtention d'un dataset ChEMBL36 préfiltré selon des critères
import sqlite3
import pandas as pd

#Connexion aux bases
conn35 = sqlite3.connect("chembl_35_sqlite/chembl_35.db")
conn36 = sqlite3.connect("chembl_36_sqlite/chembl_36.db")

#Extraire interactions ChEMBL35
pairs35 = pd.read_sql("""
SELECT DISTINCT
    td.chembl_id AS target_chembl_id,
    md.chembl_id AS molecule_chembl_id
FROM activities a
JOIN assays ON a.assay_id = assays.assay_id
JOIN target_dictionary td ON assays.tid = td.tid
JOIN molecule_dictionary md ON a.molregno = md.molregno
WHERE td.target_type = 'SINGLE PROTEIN'
""", conn35)

print("Nombre de couples uniques ChEMBL35 :", len(pairs35))
pairs35_set = set(zip(pairs35.target_chembl_id, pairs35.molecule_chembl_id))

#Extraire interactions pertinentes de ChEMBL36
query_chembl36_filtered = """
SELECT
    a.activity_id,

    -- ligand
    md.chembl_id AS molecule_chembl_id,
    cs.canonical_smiles,

    -- target
    td.chembl_id AS target_chembl_id,
    td.pref_name AS target_name,
    td.organism,
    td.target_type,

    -- sequence
    csq.accession,
    csq.sequence,

    -- activité
    a.standard_type,
    a.standard_relation,
    a.standard_value,
    a.standard_units,
    a.pchembl_value,

    -- assay
    assays.assay_id,
    assays.assay_type,
    assays.assay_organism,
    assays.description AS assay_description,

    -- publication
    docs.chembl_id AS document_chembl_id,
    docs.journal,
    docs.year,
    docs.title

FROM activities a

JOIN assays ON a.assay_id = assays.assay_id
JOIN target_dictionary td ON assays.tid = td.tid
JOIN molecule_dictionary md ON a.molregno = md.molregno

LEFT JOIN compound_structures cs
    ON md.molregno = cs.molregno

LEFT JOIN docs
    ON assays.doc_id = docs.doc_id

LEFT JOIN target_components tc
    ON td.tid = tc.tid

LEFT JOIN component_sequences csq
    ON tc.component_id = csq.component_id

WHERE td.target_type = 'SINGLE PROTEIN'
AND assays.assay_type = 'B'
AND a.standard_type IN ('Ki','Kd','IC50','EC50')
AND td.organism IN ('Homo sapiens','Escherichia coli')
AND a.pchembl_value >= 7
AND a.standard_value IS NOT NULL
"""

dataset36_filtered = pd.read_sql(query_chembl36_filtered, conn36)
print("Interactions ChEMBL36 filtrées :", dataset36_filtered.shape[0])

#Supprimer grouoes en commun avec ChEMBL35
dataset36_filtered['pair'] = list(zip(dataset36_filtered.target_chembl_id, dataset36_filtered.molecule_chembl_id))
new_dataset = dataset36_filtered[~dataset36_filtered['pair'].isin(pairs35_set)].drop(columns=['pair'])
print("Nouvelles paires ligand–protéine après comparaison :", new_dataset.shape[0])

#Sauvegarde
new_dataset.to_csv("chembl36_new_pairs_filtered_optimized.csv", index=False)



#%% Récupération des PDB associés aux Uniprot des protéines de ChEMBL36 après pré-filtrage

import sqlite3
import pandas as pd
import requests
from tqdm import tqdm

# Import de la base ChEMBL36
conn36 = sqlite3.connect("chembl_36_sqlite/chembl_36.db")

# Filtre initial
query_chembl36_filtered = """
SELECT
    a.activity_id,
    md.chembl_id AS molecule_chembl_id,
    cs.canonical_smiles,
    td.chembl_id AS target_chembl_id,
    td.pref_name AS target_name,
    td.organism,
    td.target_type,
    csq.accession,
    csq.sequence,
    a.standard_type,
    a.standard_relation,
    a.standard_value,
    a.standard_units,
    a.pchembl_value,
    assays.assay_id,
    assays.assay_type,
    assays.assay_organism,
    assays.description AS assay_description,
    docs.chembl_id AS document_chembl_id,
    docs.journal,
    docs.year,
    docs.title
FROM activities a
JOIN assays ON a.assay_id = assays.assay_id
JOIN target_dictionary td ON assays.tid = td.tid
JOIN molecule_dictionary md ON a.molregno = md.molregno
LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN docs ON assays.doc_id = docs.doc_id
LEFT JOIN target_components tc ON td.tid = tc.tid
LEFT JOIN component_sequences csq ON tc.component_id = csq.component_id
WHERE td.target_type = 'SINGLE PROTEIN'
AND assays.assay_type = 'B'
AND a.standard_type IN ('Ki','Kd','IC50','EC50')
AND td.organism IN ('Homo sapiens','Escherichia coli')
AND a.pchembl_value >= 7
AND a.standard_value IS NOT NULL
"""
dataset36_filtered = pd.read_sql(query_chembl36_filtered, conn36) # Affinité assez grande avec le score, homo sapiens, single protein, Ki/Kd/EC50/IC50

#Ajout colonne avec PDB des protéines en utilisant leur UNIPROT
import time

def fetch_pdbs_from_uniprot(uniprot_id):
    if pd.isna(uniprot_id):
        return None

    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry"
    }

    try:
        r = requests.post(url, json=query)
        if r.status_code != 200:
            return None

        data = r.json()

        if "result_set" not in data:
            return None

        pdb_ids = [item["identifier"] for item in data["result_set"]]
        return ",".join(pdb_ids) if pdb_ids else None

    except Exception as e:
        return None


uniprot_ids = dataset36_filtered['accession'].dropna().unique().tolist()

all_pdbs = {}
for uid in tqdm(uniprot_ids, desc="Récupération PDBs"):
    all_pdbs[uid] = fetch_pdbs_from_uniprot(uid)
    time.sleep(0.2)

dataset36_filtered["PDB IDs"] = dataset36_filtered["accession"].map(all_pdbs)

#Sauvegarde CSV
dataset36_filtered.to_csv("chembl_uniprot_single_pdb.csv", index=False)



#%% Application du filtre Biolip2 à Chembl36 préfiltré + enlever les Uniprot sans correspondances PDB

import pandas as pd
import json

#Chargeemnt Dataset
dataset36_filtered = pd.read_csv("ChEMBL36_filtré_PDB.csv")

#Compter nombre d'Uniprot uniques
initial_count = len(dataset36_filtered['accession'].dropna().unique())

#Supprimer les lignes sans PDB
# Créer une copie du dataset pour éviter de modifier l'original
filtered_dataset_no_pdb = dataset36_filtered[dataset36_filtered['PDB IDs'].notna()]
count_no_pdb = len(filtered_dataset_no_pdb['accession'].dropna().unique())

def extract_pdbs(obj, pdb_set):
    """Explore récursivement n'importe quelle structure JSON et extrait les PDB."""
    
    # Cas 1 : string : vérifier si c'est un PDB
    if isinstance(obj, str):
        # Un PDB est typiquement 4 caractères alphanumériques
        pdb = obj.strip().upper()  
        if len(obj) == 4 and obj.isalnum():
            pdb_set.add(obj)
    
    # Cas 2 : liste : explorer chaque élément
    elif isinstance(obj, list):
        for item in obj:
            extract_pdbs(item, pdb_set)
    
    # Cas 3 : dict : explorer chaque valeur
    elif isinstance(obj, dict):
        for key, value in obj.items():
            extract_pdbs(value, pdb_set)

    # Cas 4 : autre → ignorer
    else:
        pass

with open("PDB ID finetuning.json", "r") as file:
    data = json.load(file)

# Set pour éviter les doublons
pdb_set = set()

# Extraction récursive
extract_pdbs(data, pdb_set)
pdb_set = {p.upper() for p in pdb_set}
print(f"Nombre total de PDBs uniques trouvés : {len(pdb_set)}")

#Comparer PDB
def filter_pdb(pdb_ids):

    if pd.isna(pdb_ids):
        return False

    pdb_list = str(pdb_ids).split(",")

    for pdb in pdb_list:
        if pdb.strip().upper() in pdb_set:
            return False

    return True

# Créer un DataFrame filtré où aucun PDB du dataset n'est présent dans la liste Biolip2_article
final_filtered_dataset = filtered_dataset_no_pdb[filtered_dataset_no_pdb['PDB IDs'].apply(filter_pdb)]

# Nombre final d'UniProt IDs uniques après filtrage des PDBs présents dans la liste Biolip2_article
final_count = len(final_filtered_dataset['accession'].dropna().unique())
print(f"Nombre après élimination des lignes avec PDBs présents dans le JSON : {final_count}")

#Sauvegarde
final_filtered_dataset.to_csv("Chembl_sans_biolip2.csv", index=False)


#%% Filtre final: éviter que dataset contiennent des protéines proche de celles d'entrainement

import pandas as pd
import os
import urllib.request
import gzip
import shutil
import time
import requests
from collections import defaultdict

#Chargement des données
OLD_FILE = "PDB finetuning article.csv"
NEW_FILE = "Résultats/Résultats_finaux_propres/Biolip2_2026_vs_biolip2_article_plusieurs_ligands.csv"               # <- adapter
OUTPUT_DIR = "resultsbiolip"
IDENTITY_THRESHOLD = 0.3

SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
SIFTS_GZ = "pdb_chain_uniprot.tsv.gz"
SIFTS_TSV = "pdb_chain_uniprot.tsv"

def get_sifts():
    if not os.path.exists(SIFTS_TSV):
        if not os.path.exists(SIFTS_GZ):

            urllib.request.urlretrieve(SIFTS_URL, SIFTS_GZ)
        with gzip.open(SIFTS_GZ, "rb") as f_in, open(SIFTS_TSV, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        print("SIFTS déjà présent.")

    print("Chargement SIFTS...")
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

#Charger jeux de données
def load_old(filepath, sifts):
    df = pd.read_csv(filepath)
    df["PDB_ID"] = df["PDB_ID"].str.strip().str.lower()
    merged = df.merge(sifts, left_on="PDB_ID", right_on="PDB", how="left")
    uniprots = set(merged["SP_PRIMARY"].dropna().unique())
    return uniprots

def load_new(filepath):
    df = pd.read_csv(filepath)
    df["uniprot_id"] = df["UniProt_ID"].astype(str)
    # gérer plusieurs séparateurs
    df["uniprot_id"] = df["uniprot_id"].str.replace(",", ";")
    df["uniprot_id"] = df["uniprot_id"].str.split(";")

    df = df.explode("uniprot_id")
    df["uniprot_id"] = df["uniprot_id"].str.strip()
    uniprots = set(df["uniprot_id"].dropna().unique())
    sample_ids = sorted(uniprots)[:30]
    print(f"Total UniProt valides : {len(uniprots)}")
    print("     Sample:", list(uniprots)[:10])

    return df, uniprots

#Téléchargement séquences
import requests
import time

def download_sequences(uniprot_ids, output_fasta, label):

    ids = sorted(uniprot_ids)
    batch_size = 50
    all_fasta = []

    def fetch_batch(batch, retries=3):
        url = "https://rest.uniprot.org/uniprotkb/stream"
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        params = {"query": query, "format": "fasta"}
        for attempt in range(retries):
            try:
                r = requests.get(url, params=params, timeout=60)
                if r.ok and r.text.strip():
                    return r.text.strip()
            except Exception as e:
            time.sleep(1.5 * (attempt + 1))  # backoff exponentiel
        return None

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        result = fetch_batch(batch)
        if result:
            all_fasta.append(result)
        else:
            for acc in batch:
                single_result = fetch_batch([acc], retries=2)
                if single_result:
                    all_fasta.append(single_result)
                else:
                    print(f"Echec obtention")

        done = min(i + batch_size, len(ids))
        if done % 500 == 0 or done == len(ids):
            print(f"    {done}/{len(ids)}")

        time.sleep(2.0)

    with open(output_fasta, "w") as f:
        f.write("\n".join(all_fasta) + "\n")

    n = sum(1 for l in open(output_fasta) if l.startswith(">"))
    print(f"{n} séquences récupérées")

#Parser Fasta
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
                # >sp|P12345|NAME ... ou >P12345 ...
                header = line[1:]
                parts = header.split("|")
                if len(parts) >= 2:
                    current_id = parts[1]
                else:
                    current_id = header.split()[0]
                current_seq = []
            elif current_id:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)

    return seqs


#Calcul identitéde séquence par kmer
def kmer_identity(seq1, seq2, k=3):

    if not seq1 or not seq2:
        return 0.0

    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))

    if not kmers1 or not kmers2:
        return 0.0

    shared = len(kmers1 & kmers2)
    total = min(len(kmers1), len(kmers2))
    return shared / total


def find_homologs_kmer(old_seqs, new_seqs, threshold=0.3):
    kmer_threshold = threshold * 0.7  # conservateur

    old_list = list(old_seqs.items())
    matched = set()
    total = len(new_seqs)
    print(f"\n[4] Recherche d'homologues (k-mer, {total} queries vs {len(old_list)} targets)...")
    print(f"    Seuil k-mer = {kmer_threshold:.2f} (proxy pour {threshold*100:.0f}% identité)")

    # Pre-compute k-mers pour old
    K = 3
    old_kmers = {}
    for acc, seq in old_list:
        old_kmers[acc] = set(seq[i:i+K] for i in range(len(seq) - K + 1))

    t0 = time.time()
    for idx, (new_acc, new_seq) in enumerate(new_seqs.items()):
        if not new_seq:
            continue

        new_kset = set(new_seq[i:i+K] for i in range(len(new_seq) - K + 1))
        if not new_kset:
            continue

        n_new = len(new_kset)

        for old_acc, old_kset in old_kmers.items():
            shared = len(new_kset & old_kset)
            denom = min(n_new, len(old_kset))
            if denom > 0 and (shared / denom) >= kmer_threshold:
                matched.add(new_acc)
                break  # un seul hit suffit

        if (idx + 1) % 100 == 0 or idx + 1 == total:
            elapsed = time.time() - t0
            print(f"    {idx+1}/{total}  matched={len(matched)}  ({elapsed:.0f}s)")

    return matched

#Export
def export(truly_novel, has_homolog, novel_entries, new_df):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    p1 = os.path.join(OUTPUT_DIR, "truly_novel_uniprots.txt")
    with open(p1, "w") as f:
        for u in sorted(truly_novel):
            f.write(f"{u}\n")

    p2 = os.path.join(OUTPUT_DIR, "truly_novel_binding_sites.csv")
    novel_entries.to_csv(p2, index=False)

    p3 = os.path.join(OUTPUT_DIR, "summary.txt")
    with open(p3, "w") as f:
        f.write(f"Seuil identité (approx) : {IDENTITY_THRESHOLD*100:.0f}%\n")
        f.write(f"Total nouveau dataset   : {len(new_df)} entrées\n")
        f.write(f"UniProt nouveau         : {new_df['uniprot_id'].nunique()}\n")
        f.write(f"Avec homologue ≥30%     : {len(has_homolog)}\n")
        f.write(f"Vrais nouveaux (<30%)   : {len(truly_novel)}\n")
        f.write(f"Binding sites nouveaux  : {len(novel_entries)}\n")

    print(f"  {p1}  ({len(truly_novel)} protéines)")
    print(f"  {p2}  ({len(novel_entries)} binding sites)")
    print(f"  {p3}")



def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    old_fasta = os.path.join(OUTPUT_DIR, "old_sequences.fasta")
    new_fasta = os.path.join(OUTPUT_DIR, "new_sequences.fasta")

    # 1) SIFTS
    sifts = get_sifts()

    # 2) Charger
    old_uniprots = load_old(OLD_FILE, sifts)
    new_df, new_uniprots = load_new(NEW_FILE)

    overlap = old_uniprots & new_uniprots
    print(f"\n    Overlap direct : {len(overlap)}")
    print(f"    Pas dans ancien: {len(new_uniprots - old_uniprots)}")

    # 3) Télécharger séquences
    download_sequences(old_uniprots, old_fasta, "ancien")
    download_sequences(new_uniprots, new_fasta, "nouveau")

    # 4) Parser FASTA
    old_seqs = parse_fasta(old_fasta)
    new_seqs = parse_fasta(new_fasta)
    print(f"\n    Séquences parsées : ancien={len(old_seqs)}, nouveau={len(new_seqs)}")

    # 5) Trouver homologues (pur Python)
    matched = find_homologs_kmer(old_seqs, new_seqs, IDENTITY_THRESHOLD)

    # overlap
    uniprot_overlap = new_uniprots & old_uniprots

    # enlever overlap avant homologie
    new_unique = new_uniprots - old_uniprots
    new_seqs_filtered = {k: v for k, v in new_seqs.items() if k in new_unique}

    # homologie
    matched = find_homologs_kmer(old_seqs, new_seqs_filtered, IDENTITY_THRESHOLD)

    # split final
    homolog_only = matched
    truly_novel = new_unique - matched

    print(f"\n    RÉSULTAT :")
    print(f"    UniProt communs        : {len(uniprot_overlap)}")
    print(f"    Homologues uniquement  : {len(homolog_only)}")
    print(f"    VRAIS nouveaux         : {len(truly_novel)}")

    novel_entries = new_df[new_df["uniprot_id"].isin(truly_novel)].copy()

    # 6) Export
    df_overlap = new_df[new_df["uniprot_id"].isin(uniprot_overlap)].copy()
    df_homolog = new_df[new_df["uniprot_id"].isin(homolog_only)].copy()
    df_novel = new_df[new_df["uniprot_id"].isin(truly_novel)].copy()
    df_overlap.to_csv(os.path.join(OUTPUT_DIR, "dataset_overlap_uniprot.csv"), index=False)
    df_homolog.to_csv(os.path.join(OUTPUT_DIR, "dataset_homologs.csv"), index=False)
    df_novel.to_csv(os.path.join(OUTPUT_DIR, "dataset_truly_novel.csv"), index=False)

    print("DONE !")
    print(f"Résultats dans ./{OUTPUT_DIR}/")



if __name__ == "__main__":
    main()

#%% Créer dossiers avec homologue, nouveaux, identiques

import os
import shutil
import pandas as pd
from tqdm import tqdm

#Configuration
SOURCE_ROOT = r"/home/sixte/Documents/Helix/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip2"

CSV_HOMOLOGUES = r"Validationdatasets/results/dataset_homologs.csv"
CSV_IDENTIQUES = r"Validationdatasets/results/dataset_overlap_uniprot.csv"
CSV_NOUVEAUX   = r"Validationdatasets/results/dataset_truly_novel.csv"

OUTPUT_DIR_HOMOLOGUES = "ValidationBiolip_homologues"
OUTPUT_DIR_IDENTIQUES = "ValidationBiolip_identiques"
OUTPUT_DIR_NOUVEAUX   = "ValidationBiolip_nouveaux"

UNIPROT_COL = "UniProt_ID" 

def get_uniprot_list(csv_path, col_name):
    """Lit un CSV et retourne une liste (set) des UniProt IDs uniques."""
    if not os.path.exists(csv_path):
        return set()
    
    try:
        df = pd.read_csv(csv_path, dtype=str)
        if col_name not in df.columns:
            return set()
        # Nettoyage : enlève les espaces, met en majuscule, retire les NaN
        ids = df[col_name].dropna().str.strip().str.upper().unique()
        return set(ids)
    except Exception as e:
        return set()

def organize_folders(uniprot_ids, source_root, dest_root, category_name):
    if not uniprot_ids:
        return
    # Création du dossier de destination
    os.makedirs(dest_root, exist_ok=True)
    # Liste des dossiers existants dans la source
    existing_source_dirs = set(os.listdir(source_root))
    
    copied_count = 0
    missing_count = 0
        
    for uniprot_id in tqdm(uniprot_ids, desc=f"Copie {category_name}", unit="dossier"):
        src_path = os.path.join(source_root, uniprot_id)
        dst_path = os.path.join(dest_root, uniprot_id)
        
        if uniprot_id not in existing_source_dirs or not os.path.isdir(src_path):
            missing_count += 1
            # print(f"   ⚠️  Dossier source introuvable pour : {uniprot_id}")
            continue
        
        # Si la destination existe déjà, on la saute (ou on pourrait l'écraser avec shutil.rmtree avant)
        if os.path.exists(dst_path):
            continue
            
        try:
            shutil.copytree(src_path, dst_path)
            copied_count += 1
        except Exception as e:
            print(f"Erreur copie de {uniprot_id}: {e}")

def main():    
    if not os.path.exists(SOURCE_ROOT):
        print(f"Le dossier source n'existe pas : {SOURCE_ROOT}")
        return

    # 1. Extraction des listes d'UniProt
    set_homologues = get_uniprot_list(CSV_HOMOLOGUES, UNIPROT_COL)
    set_identiques = get_uniprot_list(CSV_IDENTIQUES, UNIPROT_COL)
    set_nouveaux   = get_uniprot_list(CSV_NOUVEAUX,   UNIPROT_COL)
    
    print(f"   - Homologues : {len(set_homologues)} IDs")
    print(f"   - Identiques : {len(set_identiques)} IDs")
    print(f"   - Nouveaux   : {len(set_nouveaux)} IDs")

    # 2. Organisation des dossiers    
    organize_folders(set_homologues, SOURCE_ROOT, OUTPUT_DIR_HOMOLOGUES, "HOMOLOGUES")
    organize_folders(set_identiques, SOURCE_ROOT, OUTPUT_DIR_IDENTIQUES, "IDENTIQUES")
    organize_folders(set_nouveaux,   SOURCE_ROOT, OUTPUT_DIR_NOUVEAUX,   "NOUVEAUX")

if __name__ == "__main__":
    main()


#%% Resource deadlock avoided

import os
import shutil
import pandas as pd
from tqdm import tqdm
import stat

# ============================================================
#                   CONFIGURATION
# ============================================================

SOURCE_ROOT = r"/home/sixte/Documents/Helix/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationEcoli_sansKiKd"
CSV_HOMOLOGUES = r"Validationdatasets/results/dataset_homologs.csv"
CSV_IDENTIQUES = r"Validationdatasets/results/dataset_overlap_uniprot.csv"
CSV_NOUVEAUX   = r"Validationdatasets/results/dataset_truly_novel.csv"

OUTPUT_DIR_HOMOLOGUES = "ValidationBiolip_homologues"
OUTPUT_DIR_IDENTIQUES = "ValidationBiolip_identiques"
OUTPUT_DIR_NOUVEAUX   = "ValidationBiolip_nouveaux"

UNIPROT_COL = "UniProt_ID" 

# On force la réécriture pour être sûr d'avoir des dossiers propres
FORCE_OVERWRITE = True 

# ============================================================
#                   FONCTION DE COPIE "BUNKER" (Ultra-Sécurisée)
# ============================================================

def safe_copy_tree(src, dst):    
    # 1. Création du dossier de destination
    try:
        if not os.path.exists(dst):
            os.makedirs(dst)
    except OSError as e:
        print(f"Impossible de créer le dossier {dst}: {e}")
        return

    # 2. Liste du contenu source (sécurisée)
    try:
        items = os.listdir(src)
    except OSError as e:
        # Si on ne peut même pas lister le dossier source (Deadlock ici), on abandonne ce dossier
        print(f"Impossible de lire le dossier source {src} (Deadlock/Permission). Dossier ignoré.")
        return
    except PermissionError:
        print(f"Permission refusée pour lire : {src}")
        return

    for item in items:
        src_path = os.path.join(src, item)
        dst_path = os.path.join(dst, item)
        # islink renvoie True si c'est un lien, même s'il pointe vers un dossier.
        if os.path.islink(src_path):
            continue

        #Traitement des dossiers
        if os.path.isdir(src_path):
            safe_copy_tree(src_path, dst_path)
            
        #Traitement des fichiers
        elif os.path.isfile(src_path):
            try:
                shutil.copy2(src_path, dst_path)
            except OSError as e:
                if e.errno == 35:
                    # Deadlock sur un fichier spécifique -> on ignore
                    pass 
                else:
                    # Autre erreur (espace disque, permission) -> on ignore aussi pour continuer
                    pass
            except Exception:
                pass

#Lecture
def get_uniprot_list(csv_path, col_name):
    if not os.path.exists(csv_path):
        return set()
    try:
        df = pd.read_csv(csv_path, dtype=str)
        if col_name not in df.columns:
            return set()
        ids = df[col_name].dropna().str.strip().str.upper().unique()
        return set(ids)
    except Exception as e:
        return set()

def organize_folders(uniprot_ids, source_root, dest_root, category_name):
    if not uniprot_ids:
        return

    # Création du dossier racine de destination avec gestion d'erreur stricte
    try:
        os.makedirs(dest_root, exist_ok=True)
    except OSError as e:
        print(f"Erreur création dossier racine")
        return

    try:
        existing_source_dirs = set(os.listdir(source_root))
    except FileNotFoundError:
        return
    except OSError as e:
        return
    
    copied_count = 0
    missing_count = 0
    error_count = 0
    skipped_count = 0
    
    if FORCE_OVERWRITE:
    
    for uniprot_id in tqdm(uniprot_ids, desc=f"Copie {category_name}", unit="dossier"):
        src_path = os.path.join(source_root, uniprot_id)
        dst_path = os.path.join(dest_root, uniprot_id)
        
        if uniprot_id not in existing_source_dirs or not os.path.isdir(src_path):
            missing_count += 1
            continue
        
        if os.path.exists(dst_path) or os.path.islink(dst_path):
            if FORCE_OVERWRITE:
                try:
                    shutil.rmtree(dst_path)
                except Exception as e:
                    error_count += 1
                    continue
            else:
                skipped_count += 1
                continue
            
        try:
            safe_copy_tree(src_path, dst_path)
            copied_count += 1
        except Exception as e:
            error_count += 1

def main():
    
    if not os.path.exists(SOURCE_ROOT):
        return

    set_homologues = get_uniprot_list(CSV_HOMOLOGUES, UNIPROT_COL)
    set_identiques = get_uniprot_list(CSV_IDENTIQUES, UNIPROT_COL)
    set_nouveaux   = get_uniprot_list(CSV_NOUVEAUX,   UNIPROT_COL)
    
    print(f"   - Homologues : {len(set_homologues)} IDs")
    print(f"   - Identiques : {len(set_identiques)} IDs")
    print(f"   - Nouveaux   : {len(set_nouveaux)} IDs")

    # Lancement séquentiel
    organize_folders(set_homologues, SOURCE_ROOT, OUTPUT_DIR_HOMOLOGUES, "HOMOLOGUES")
    organize_folders(set_identiques, SOURCE_ROOT, OUTPUT_DIR_IDENTIQUES, "IDENTIQUES")
    organize_folders(set_nouveaux,   SOURCE_ROOT, OUTPUT_DIR_NOUVEAUX,   "NOUVEAUX")

if __name__ == "__main__":
    main()