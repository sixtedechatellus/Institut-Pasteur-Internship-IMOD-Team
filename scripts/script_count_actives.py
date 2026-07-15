#!/usr/bin/env python3
"""
Calcule la moyenne et l'écart-type du nombre d'actifs par protéine
dans tous les sous-dossiers contenant un fichier decoys_final.mol2.

Usage :
    python compute_mean_actives.py <dossier_racine>
"""

import os
import sys
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python compute_mean_actives.py <dossier_racine>")
    sys.exit(1)

root = sys.argv[1]
if not os.path.exists(root):
    print(f"❌ Dossier inexistant : {root}")
    sys.exit(1)

counts = []   # liste du nombre d'actifs pour chaque dossier traité

for subdir, _, files in os.walk(root):
    if "actives_final.ism" in files:
        actives_path = os.path.join(subdir, "actives_final.ism")
        if os.path.exists(actives_path):
            with open(actives_path) as fin:
                n_actives = sum(1 for line in fin if line.strip())
            counts.append(n_actives)
            print(f"[✔] {os.path.basename(subdir)} : {n_actives} actifs")
        else:
            print(f"[⚠] {os.path.basename(subdir)} : pas de actives_final.ism")

if not counts:
    print("Aucun dossier avec decoys_final.mol2 trouvé.")
    sys.exit(0)

counts = np.array(counts)
mean_actives = counts.mean()
std_actives = counts.std(ddof=1)  # écart-type échantillon

print("\n================ Résumé ================")
print(f"Nombre de dossiers avec decoys_final.mol2 : {len(counts)}")
print(f"Total d'actifs comptés : {counts.sum()}")
print(f"Moyenne d'actifs par protéine : {mean_actives:.2f}")
print(f"Écart-type du nombre d'actifs : {std_actives:.2f}")
