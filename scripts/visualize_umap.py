"""
visualize_umap.py — avec étiquetage numéroté des ligands actifs
Usage : python visualize_umap.py --emb_dir ./resultats/embeddings
                                  --scores_dir ./dossier_scores
                                  --targets prot1 prot2
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import umap

# ─────────────────────────────────────────────
# Paramètres
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--emb_dir",    type=str, default="./resultats/embeddings")
parser.add_argument("--scores_dir", type=str, default=None)
parser.add_argument("--targets",    type=str, nargs="+", default=None)
parser.add_argument("--output",     type=str, default="umap_latent_space.png")
args = parser.parse_args()

# ─────────────────────────────────────────────
# 1. Chargement
# ─────────────────────────────────────────────
emb_root = args.emb_dir

if args.targets is None:
    targets = sorted([
        d for d in os.listdir(emb_root)
        if os.path.isdir(os.path.join(emb_root, d))
    ])
else:
    targets = args.targets

print(f"[INFO] Proteines : {targets}")

cmap      = cm.get_cmap("tab10", max(len(targets), 1))
color_map = {prot: cmap(i) for i, prot in enumerate(targets)}

all_embeddings = []
all_types      = []
all_proteins   = []
all_smiles     = []
active_labels  = []    # [{num, protein, smiles}, ...]

# ─────────────────────────────────────────────
# Fonction de chargement SMILES robuste
# ─────────────────────────────────────────────
def charger_smiles(scores_dir, protein, mol_reps, labels):
    """
    Retourne un tableau smiles_array de longueur len(mol_reps).
    Stratégie :
      1. Alignement exact (même longueur CSV et embeddings)
      2. Alignement par les actifs seulement (si longueurs diffèrent)
      3. Fallback : '?' pour tous
    """
    smiles_array = np.array(["?"] * len(mol_reps), dtype=object)

    if scores_dir is None:
        print(f"    [WARN] --scores_dir non fourni → SMILES = '?'")
        return smiles_array

    # ── Trouver le CSV correspondant ──────────────────────────
    try:
        all_files    = os.listdir(scores_dir)
    except FileNotFoundError:
        print(f"    [WARN] Dossier scores introuvable : {scores_dir}")
        return smiles_array

    csv_candidates = [f for f in all_files if f.endswith(".csv") and protein in f]

    if not csv_candidates:
        print(f"    [WARN] Aucun CSV trouvé pour '{protein}' dans {scores_dir}")
        return smiles_array

    csv_path = os.path.join(scores_dir, csv_candidates[0])
    print(f"    → CSV trouvé : {csv_candidates[0]}")

    try:
        df = pd.read_csv(csv_path, engine='python', on_bad_lines='skip')
        df.columns = df.columns.str.strip().str.lower()

        # ── Debug colonnes ─────────────────────────────────────
        print(f"    → Colonnes CSV     : {df.columns.tolist()}")
        print(f"    → Lignes CSV       : {len(df)}")
        print(f"    → Embeddings       : {len(mol_reps)}")

        if 'ligand' not in df.columns:
            print(f"    [WARN] Colonne 'ligand' absente !")
            print(f"           Colonnes disponibles : {df.columns.tolist()}")
            return smiles_array

        # ── Stratégie 1 : alignement exact ────────────────────
        if len(df) == len(mol_reps):
            smiles_array = df['ligand'].values.astype(str)
            print(f"    → SMILES chargés ✓  (alignement exact, {len(df)} lignes)")
            return smiles_array

        # ── Stratégie 2 : alignement actifs uniquement ─────────
        print(f"    [WARN] Longueurs différentes ({len(df)} CSV ≠ {len(mol_reps)} emb)")
        print(f"           Tentative d'alignement par les actifs…")

        df_clean = df.dropna(subset=['label']).copy()
        df_clean['label'] = df_clean['label'].astype(int)

        actifs_csv    = df_clean[df_clean['label'] == 1]['ligand'].values.astype(str)
        n_actifs_emb  = int((labels == 1).sum())

        print(f"    → Actifs dans CSV  : {len(actifs_csv)}")
        print(f"    → Actifs emb label : {n_actifs_emb}")

        if len(actifs_csv) == n_actifs_emb:
            # On place les SMILES aux positions des actifs dans mol_reps
            active_positions = np.where(labels == 1)[0]
            for pos, smi in zip(active_positions, actifs_csv):
                smiles_array[pos] = smi
            print(f"    → SMILES actifs chargés ✓ "
                  f"(alignement par label, {len(actifs_csv)} actifs)")
        else:
            print(f"    [WARN] Impossible d'aligner les actifs "
                  f"({len(actifs_csv)} CSV ≠ {n_actifs_emb} emb) → SMILES = '?'")

    except Exception as e:
        print(f"    [WARN] Erreur lecture CSV — {e}")

    return smiles_array


# ─────────────────────────────────────────────
# Boucle principale
# ─────────────────────────────────────────────
for protein in targets:
    folder = os.path.join(emb_root, protein)

    needed = ["mol_reps.npy", "pocket_reps.npy", "labels.npy"]
    if not all(os.path.exists(os.path.join(folder, f)) for f in needed):
        print(f"[WARN] Fichiers manquants pour {protein} — ignoré.")
        continue

    mol_reps    = np.load(os.path.join(folder, "mol_reps.npy"))
    pocket_reps = np.load(os.path.join(folder, "pocket_reps.npy"))
    labels      = np.load(os.path.join(folder, "labels.npy"))

    print(f"\n── {protein} ──────────────────────────────────────")

    # ── SMILES ──────────────────────────────────────────────
    smiles_array = charger_smiles(args.scores_dir, protein, mol_reps, labels)

    # ── Poche ────────────────────────────────────────────────
    pocket_mean = pocket_reps.mean(axis=0, keepdims=True)
    all_embeddings.append(pocket_mean)
    all_types.append("pocket")
    all_proteins.append(protein)
    all_smiles.append("")

    # ── Vrais ligands ────────────────────────────────────────
    active_mask   = (labels == 1)
    actifs        = mol_reps[active_mask]
    actifs_smiles = smiles_array[active_mask]

    if len(actifs) > 0:
        all_embeddings.append(actifs)
        all_types.extend(["active"] * len(actifs))
        all_proteins.extend([protein]  * len(actifs))

        for s in actifs_smiles:
            num = len(active_labels) + 1
            all_smiles.append(str(s))
            active_labels.append({
                "num"    : num,
                "protein": protein,
                "smiles" : str(s),
            })

    # ── Decoys ───────────────────────────────────────────────
    decoys = mol_reps[labels == 0]
    if len(decoys) > 0:
        all_embeddings.append(decoys)
        all_types.extend(["decoy"]   * len(decoys))
        all_proteins.extend([protein] * len(decoys))
        all_smiles.extend([""]        * len(decoys))

    print(f"    Résumé : {len(actifs)} actifs, {len(decoys)} decoys")
    for i, entry in enumerate([e for e in active_labels if e["protein"] == protein]):
        smi_preview = entry["smiles"][:60] + "…" if len(entry["smiles"]) > 60 else entry["smiles"]
        print(f"      Actif #{entry['num']} : {smi_preview}")

# ── Conversion finale ────────────────────────────────────────
all_embeddings = np.concatenate(all_embeddings, axis=0)
all_types      = np.array(all_types,    dtype=str)
all_proteins   = np.array(all_proteins, dtype=str)
all_smiles     = np.array(all_smiles,   dtype=str)

active_idx = np.where(all_types == "active")[0]

print(f"\n[INFO] Total embeddings : {len(all_embeddings)}")
print(f"       Poches  : {(all_types == 'pocket').sum()}")
print(f"       Actifs  : {(all_types == 'active').sum()}")
print(f"       Decoys  : {(all_types == 'decoy').sum()}")
print(f"       SMILES '?' restants : {(all_smiles == '?').sum()}")

# ─────────────────────────────────────────────
# 2. UMAP
# ─────────────────────────────────────────────
print("\n[INFO] Calcul UMAP…")
reducer = umap.UMAP(
    n_components=2,
    n_neighbors=15,
    min_dist=0.1,
    metric="cosine",
    random_state=42
)
emb_2d = reducer.fit_transform(all_embeddings)
print("[INFO] UMAP terminé.")

# ─────────────────────────────────────────────
# 3. Visualisation
# ─────────────────────────────────────────────
marker_map = {"pocket": "*", "active": "o", "decoy": "x"}
size_map   = {"pocket": 400, "active": 80,  "decoy": 25}
alpha_map  = {"pocket": 1.0, "active": 0.85,"decoy": 0.35}
zorder_map = {"pocket": 3,   "active": 2,   "decoy": 1}

n_actifs   = len(active_labels)
show_table = (n_actifs > 0) and (n_actifs <= 35)

if show_table:
    fig      = plt.figure(figsize=(14, 13))
    ax       = fig.add_axes([0.05, 0.32, 0.90, 0.64])
    ax_table = fig.add_axes([0.02, 0.01, 0.96, 0.28])
else:
    fig, ax  = plt.subplots(figsize=(12, 9))
    ax_table = None

# ── Scatter ──────────────────────────────────────────────────
for mol_type in ["decoy", "active", "pocket"]:
    mask = (all_types == mol_type)
    if mask.sum() == 0:
        continue

    xs     = emb_2d[mask, 0]
    ys     = emb_2d[mask, 1]
    colors = [color_map[p] for p in all_proteins[mask]]

    ax.scatter(
        xs, ys,
        c=colors,
        marker=marker_map[mol_type],
        s=size_map[mol_type],
        alpha=alpha_map[mol_type],
        edgecolors="black" if mol_type == "pocket" else "none",
        linewidths=1.5,
        zorder=zorder_map[mol_type],
    )

# ── Numérotation des actifs ──────────────────────────────────
for i, idx in enumerate(active_idx):
    num  = i + 1
    x, y = emb_2d[idx, 0], emb_2d[idx, 1]
    prot = all_proteins[idx]

    ax.annotate(
        str(num),
        xy=(x, y),
        fontsize=7, fontweight="bold",
        color="black",
        ha="center", va="center",
        zorder=6,
        bbox=dict(
            boxstyle="circle,pad=0.18",
            facecolor="white",
            edgecolor=color_map[prot],
            linewidth=1.5,
            alpha=0.90,
        )
    )

# ── Légendes ─────────────────────────────────────────────────
legend_types = [
    plt.Line2D([0],[0], marker="*", color="gray", markersize=14,
               linestyle="None", markeredgecolor="black",
               label="Poche protéique"),
    plt.Line2D([0],[0], marker="o", color="gray", markersize=8,
               linestyle="None", label="Vrai ligand (actif, numéroté)"),
    plt.Line2D([0],[0], marker="x", color="gray", markersize=7,
               linestyle="None", label="Decoy"),
]
legend_proteins = [
    mpatches.Patch(color=color_map[prot], label=prot)
    for prot in targets
]

leg1 = ax.legend(handles=legend_types,
                 loc="upper left", fontsize=9,
                 title="Type", framealpha=0.9)
ax.add_artist(leg1)
ax.legend(handles=legend_proteins,
          loc="upper right", fontsize=8,
          title="Protéine", framealpha=0.9,
          ncol=max(1, len(targets) // 15))

ax.set_title(
    "Espace latent DrugCLIP\n"
    "Poches, ligands actifs (⊙ numérotés) et decoys (UMAP cosine)",
    fontsize=13, fontweight="bold"
)
ax.set_xlabel("UMAP 1", fontsize=11)
ax.set_ylabel("UMAP 2", fontsize=11)

# ── Table des ligands ────────────────────────────────────────
if ax_table is not None and active_labels:
    ax_table.axis("off")

    table_data = []
    for entry in active_labels:
        smiles_trunc = (entry["smiles"][:90] + "…"
                        if len(entry["smiles"]) > 90
                        else entry["smiles"])
        table_data.append([str(entry["num"]), entry["protein"], smiles_trunc])

    tbl = ax_table.table(
        cellText  = table_data,
        colLabels = ["#", "Protéine", "SMILES (tronqué à 90 car.)"],
        cellLoc   = "left",
        loc       = "upper center",
        colWidths = [0.035, 0.10, 0.865],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6.5)
    tbl.scale(1, 1.25)

    for col in range(3):
        tbl[0, col].set_facecolor("#444444")
        tbl[0, col].set_text_props(color="white", fontweight="bold")

    for row_idx, entry in enumerate(active_labels, start=1):
        base_color = color_map[entry["protein"]]
        cell_color = (*base_color[:3], 0.18)
        for col in range(3):
            tbl[row_idx, col].set_facecolor(cell_color)

    ax_table.set_title(
        "Correspondance ligands actifs numérotés → SMILES",
        fontsize=9, fontweight="bold", pad=3, loc="left"
    )

elif active_labels:
    ax.text(0.01, 0.01,
            f"⚠ {n_actifs} ligands actifs numérotés\n"
            f"→ voir {args.output.replace('.png', '_ligands.csv')}",
            transform=ax.transAxes, fontsize=8,
            va="bottom", ha="left",
            bbox=dict(facecolor="lightyellow", edgecolor="gray",
                      alpha=0.9, boxstyle="round,pad=0.4"))

# ─────────────────────────────────────────────
# 4. Sauvegarde
# ─────────────────────────────────────────────
plt.savefig(args.output, dpi=300, bbox_inches="tight")
print(f"\n[INFO] Figure sauvegardée : {args.output}")

if active_labels:
    df_map  = pd.DataFrame(active_labels)[["num", "protein", "smiles"]]
    csv_out = args.output.replace(".png", "_ligands.csv")
    df_map.to_csv(csv_out, index=False)
    print(f"[INFO] Mapping ligands  : {csv_out}")
    print("\n" + df_map.to_string(index=False))

plt.show()
