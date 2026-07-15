import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.metrics import roc_curve, auc
from pathlib import Path

DATA_DIR = Path("CSV_scores/")

all_dfs = []

for filepath in sorted(DATA_DIR.glob("*.csv")):
    prot_name = filepath.stem

    try:
        df = pd.read_csv(filepath, engine='python', on_bad_lines='skip')
        df.columns = df.columns.str.strip().str.lower()

        missing = [c for c in ['ligand','score','label'] if c not in df.columns]
        if missing:
            print(f"⚠️  {prot_name} — colonnes manquantes : {missing}")
            continue

        # ── Nettoyage NaN ──────────────────────────────────────
        n_avant = len(df)
        df = df.dropna(subset=['label', 'score'])          # supprime les NaN
        df['label'] = df['label'].astype(int)              # force en entier
        df['score'] = df['score'].astype(float)            # force en float
        n_apres = len(df)

        if n_avant != n_apres:
            print(f"⚠️  {prot_name} : {n_avant - n_apres} lignes NaN supprimées")

        df["protein"] = prot_name
        all_dfs.append(df)

        print(f"  {prot_name:35s}  "
              f"{int(df['label'].sum())} ligands / "
              f"{int((df['label']==0).sum())} leurres")

    except Exception as e:
        print(f"❌ {prot_name} — ERREUR : {e}")
        continue

all_data = pd.concat(all_dfs, ignore_index=True)
all_data.to_csv("all_scores_concatenated.csv", index=False)
print(f"\n✅ Total : {len(all_data)} entrées, "
      f"{all_data['protein'].nunique()} protéines\n")

# ============================================================
# CALCUL DES COURBES ROC PAR PROTÉINE
# ============================================================
results  = {}
ignorees = []

for prot, group in all_data.groupby("protein"):
    y_true  = group["label"].values.astype(int)
    y_score = group["score"].values.astype(float)

    # Vérifications
    if len(np.unique(y_true)) < 2:
        print(f"⚠️  {prot} : une seule classe présente, ignorée")
        ignorees.append(prot)
        continue
    if np.isnan(y_score).any() or np.isnan(y_true).any():
        print(f"⚠️  {prot} : NaN résiduels, ignorée")
        ignorees.append(prot)
        continue

    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc     = auc(fpr, tpr)

    results[prot] = {
        "fpr"  : fpr,
        "tpr"  : tpr,
        "auc"  : roc_auc,
        "n_pos": int(y_true.sum()),
        "n_neg": int((y_true == 0).sum()),
    }

print(f"\n✅ {len(results)} protéines tracées, "
      f"{len(ignorees)} ignorées ({len(ignorees)} cibles à 1 seule classe)")

results = dict(sorted(results.items(),
                      key=lambda x: x[1]["auc"], reverse=True))

# ============================================================
# FIGURE : Courbes ROC
# ============================================================
n    = len(results)
cmap = cm.get_cmap("tab20", max(n, 1))

fig, ax = plt.subplots(figsize=(10, 7))

for i, (prot, v) in enumerate(results.items()):
    ax.plot(v["fpr"], v["tpr"],
            color=cmap(i), lw=1.8,
            label=f"{prot}  (AUC = {v['auc']:.3f})")

ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="Aléatoire")

ax.set_xlabel("Taux de faux positifs (FPR)", fontsize=13)
ax.set_ylabel("Taux de vrais positifs (TPR)", fontsize=13)
ax.set_title(f"Courbes ROC — DrugCLIP ({n} protéines)", fontsize=14)
ax.legend(fontsize=7.5, loc="lower right",
          framealpha=0.85, ncol=1 if n <= 12 else 2)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.02)
ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig("roc_curves_drugclip.png", dpi=150)
plt.show()

# ============================================================
# RÉSUMÉ
# ============================================================
summary = pd.DataFrame([
    {"proteine": p, "auc": v["auc"],
     "n_ligands": v["n_pos"], "n_leurres": v["n_neg"]}
    for p, v in results.items()
])
print("\n", summary.to_string(index=False))
summary.to_csv("auc_summary.csv", index=False)
#%% autre représentation


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.metrics import roc_curve, auc
from pathlib import Path

# ============================================================
# PARAMÈTRES
# ============================================================
DATA_DIR = Path("CSV_scores/")

# ============================================================
# LECTURE + NETTOYAGE
# ============================================================
all_dfs = []

for filepath in sorted(DATA_DIR.glob("*.csv")):
    prot_name = filepath.stem

    try:
        df = pd.read_csv(filepath, engine='python', on_bad_lines='skip')
        df.columns = df.columns.str.strip().str.lower()

        missing = [c for c in ['ligand','score','label'] if c not in df.columns]
        if missing:
            print(f"⚠️  {prot_name} — colonnes manquantes : {missing}")
            continue

        n_avant = len(df)
        df = df.dropna(subset=['label', 'score'])
        df['label'] = df['label'].astype(int)
        df['score'] = df['score'].astype(float)
        n_apres = len(df)

        if n_avant != n_apres:
            print(f"⚠️  {prot_name} : {n_avant - n_apres} lignes NaN supprimées")

        df["protein"] = prot_name
        all_dfs.append(df)

        print(f"  {prot_name:35s}  "
              f"{int(df['label'].sum())} ligands / "
              f"{int((df['label']==0).sum())} leurres")

    except Exception as e:
        print(f"❌ {prot_name} — ERREUR : {e}")
        continue

all_data = pd.concat(all_dfs, ignore_index=True)
all_data.to_csv("all_scores_concatenated.csv", index=False)
print(f"\n✅ Total : {len(all_data)} entrées, "
      f"{all_data['protein'].nunique()} protéines\n")

# ============================================================
# CALCUL DES COURBES ROC PAR PROTÉINE
# ============================================================
results  = {}
ignorees = []

for prot, group in all_data.groupby("protein"):
    y_true  = group["label"].values.astype(int)
    y_score = group["score"].values.astype(float)

    if len(np.unique(y_true)) < 2:
        ignorees.append(prot)
        continue
    if np.isnan(y_score).any() or np.isnan(y_true).any():
        ignorees.append(prot)
        continue

    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc     = auc(fpr, tpr)

    results[prot] = {
        "fpr"  : fpr,
        "tpr"  : tpr,
        "auc"  : roc_auc,
        "n_pos": int(y_true.sum()),
        "n_neg": int((y_true == 0).sum()),
    }

print(f"✅ {len(results)} protéines calculées, "
      f"{len(ignorees)} ignorées (une seule classe)\n")

results = dict(sorted(results.items(),
                      key=lambda x: x[1]["auc"], reverse=True))

# ============================================================
# INTERPOLATION SUR GRILLE COMMUNE
# ============================================================
mean_fpr = np.linspace(0, 1, 500)
all_tprs = []

for prot, v in results.items():
    interp_tpr = np.interp(mean_fpr, v["fpr"], v["tpr"])
    all_tprs.append(interp_tpr)

all_tprs = np.array(all_tprs)

mean_tpr = all_tprs.mean(axis=0)
std_tpr  = all_tprs.std(axis=0)
mean_auc = np.mean([v["auc"] for v in results.values()])
std_auc  = np.std( [v["auc"] for v in results.values()])
median_auc = np.median([v["auc"] for v in results.values()])

# ============================================================
# FIGURE : courbes individuelles + moyenne + enveloppe
# ============================================================
fig, ax = plt.subplots(figsize=(8, 7))

# Courbes individuelles en gris transparent
for tpr_interp in all_tprs:
    ax.plot(mean_fpr, tpr_interp,
            color="gray", lw=0.6, alpha=0.12)

# Enveloppe ± 1 std
ax.fill_between(mean_fpr,
                np.clip(mean_tpr - std_tpr, 0, 1),
                np.clip(mean_tpr + std_tpr, 0, 1),
                color="steelblue", alpha=0.2,
                label="± 1 std")

# Courbe moyenne
ax.plot(mean_fpr, mean_tpr,
        color="steelblue", lw=2.5,
        label=f"Moyenne  (AUC = {mean_auc:.3f} ± {std_auc:.3f})")

# Ligne aléatoire
ax.plot([0, 1], [0, 1],
        "k--", lw=1, alpha=0.5, label="Aléatoire (AUC = 0.500)")

# Annotations
ax.text(0.98, 0.08,
        f"n = {len(results)} protéines\n"
        f"Médiane AUC = {median_auc:.3f}\n"
        f"{len(ignorees)} protéines ignorées",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="lightgray", alpha=0.8))

ax.set_xlabel("Taux de faux positifs (FPR)", fontsize=13)
ax.set_ylabel("Taux de vrais positifs (TPR)", fontsize=13)
ax.set_title(f"Courbes ROC — DrugCLIP\n"
             f"({len(results)} protéines DUDE, courbes individuelles en gris)",
             fontsize=13)
ax.legend(fontsize=11, loc="lower right", framealpha=0.85)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.02)
ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig("roc_mean_envelope.png", dpi=150)
plt.show()

# ============================================================
# RÉSUMÉ CSV
# ============================================================
summary = pd.DataFrame([
    {"proteine" : p,
     "auc"      : v["auc"],
     "n_ligands": v["n_pos"],
     "n_leurres": v["n_neg"]}
    for p, v in results.items()
]).sort_values("auc", ascending=False)

print(summary.to_string(index=False))
summary.to_csv("auc_summary.csv", index=False)
print(f"\nFichiers sauvegardés : roc_mean_envelope.png, auc_summary.csv")