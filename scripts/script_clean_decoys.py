#!/usr/bin/env python3
import os
import sys

MIN_SMILES_LEN = 4

if len(sys.argv) != 2:
    print("Use: python remove_identical_decoys.py <dossier_contenant_les_decoys>")
    sys.exit(1)

root = sys.argv[1]
if not os.path.exists(root):
    sys.exit(1)

for subdir, _, files in os.walk(root):
    for f in files:
        if f == "decoys_final.ism":
            path = os.path.join(subdir, f)

            cleaned = []
            removed_identical = 0
            removed_duplicates = 0
            removed_small = 0
            seen_decoys = set()

            with open(path) as fin:
                for line in fin:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) < 2:
                        cleaned.append(line + "\n")
                        continue

                    active, decoy = parts[0].strip(), parts[1].strip()

                    # Case 1 : active == decoy
                    if active == decoy:
                        removed_identical += 1
                        continue

                    # Case 2 : known decoy
                    if decoy in seen_decoys:
                        removed_duplicates += 1
                        continue
                    
                    # Case 3 : too small
                    if len(decoy) < MIN_SMILES_LEN:
                        removed_small += 1
                        continue

                    # We keep the line
                    seen_decoys.add(decoy)
                    cleaned.append(line + "\n")

            with open(path, "w") as fout:
                fout.writelines(cleaned)

            print(f"   ✅ {len(cleaned)} kept | {removed_identical} doubles removed | {removed_small} small molecules removed | {removed_duplicates} doublons supprimés")
