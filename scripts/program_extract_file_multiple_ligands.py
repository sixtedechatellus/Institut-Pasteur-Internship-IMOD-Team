import os
import shutil

src_root = "data/ValidationBiolip2"
dst_root = "data/ValidationBiolip3"

# créer le dossier de sortie
os.makedirs(dst_root, exist_ok=True)

kept = 0
skipped = 0

for subdir in os.listdir(src_root):
    subpath = os.path.join(src_root, subdir)

    if not os.path.isdir(subpath):
        continue

    keep = False

    #cas 1 : dossier ligand
    ligand_path = os.path.join(subpath, "ligands")
    if os.path.isdir(ligand_path):
        files = [f for f in os.listdir(ligand_path) if not f.startswith(".")]
        if len(files) > 1:
            keep = True

    #cas 2 : fichier actives_smiles.ism
    ism_path = os.path.join(subpath, "actives_final.ism")
    if os.path.exists(ism_path):
        with open(ism_path, "r") as f:
            lines = [l for l in f if l.strip()]
        if len(lines) > 1:
            keep = True

    #copie
    if keep:
        dst_path = os.path.join(dst_root, subdir)
        print(f" KEEP: {subdir}")
        shutil.copytree(subpath, dst_path, dirs_exist_ok=True)
        kept += 1
    else:
        print(f" SKIP: {subdir}")
        skipped += 1

print(f"Dossiers gardés: {kept}")
print(f"Dossiers ignorés: {skipped}")
