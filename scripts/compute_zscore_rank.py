import os
import sys
import pandas as pd
import numpy as np

# -----------------------------
# ARGUMENT CHECK
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python compute_zscore_rank.py <dossier_contenant_les_CSV>")
    sys.exit(1)

root = sys.argv[1]

if not os.path.exists(root):
    print(f"Erreur : dossier {root} inexistant")
    sys.exit(1)


# -----------------------------
# MAIN LOOP
# -----------------------------
rows = []

for fname in os.listdir(root):

    if not fname.endswith(".csv"):
        continue

    path = os.path.join(root, fname)

    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"⚠️ Impossible de lire {fname}: {e}")
        continue

    # nettoyage colonnes
    df.columns = df.columns.str.strip().str.lower()

    # vérification minimale
    required = {"label", "score", "ligand"}
    if not required.issubset(df.columns):
        print(f"⚠️ Colonnes manquantes dans {fname}: {df.columns.tolist()}")
        continue

    # protéine
    prot = df["protein"].iloc[0] if "protein" in df.columns else fname.replace(".csv", "")

    # tri global (important pour le rank)
    df_sorted = df.sort_values("score", ascending=False).reset_index(drop=True)

    # rangs (1-based)
    df_sorted["rank"] = df_sorted.index + 1

    # garder uniquement les actifs
    actives = df_sorted[df_sorted["label"] == 1]

    # stocker chaque actif = une ligne
    for _, row in actives.iterrows():
        rows.append({
            "PDB_ID": prot,
            "SMILES": row["ligand"],
            "Score": row["score"],   # ton "zscore"
            "Rank": row["rank"]
        })


# -----------------------------
# SAVE OUTPUT
# -----------------------------
res = pd.DataFrame(rows)

out_path = os.path.join(root, "actives_detailed.csv")
res.to_csv(out_path, index=False)

print(f"\n✅ Tableau détaillé enregistré dans {out_path}")
print(res.head(10))
