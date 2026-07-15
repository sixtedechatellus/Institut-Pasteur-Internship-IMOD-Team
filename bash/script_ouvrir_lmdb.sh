#!/bin/bash
#SBATCH --job-name=lecture_lmdb
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o logs/lecture_lmdb_%j.out
#SBATCH -e logs/lecture_lmdb_%j.err

module purge
module load cuda/11.8

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/

python3 - << "PY"
import os, lmdb, pickle

path = "data/ValidationBiolip/1A4M/mols.lmdb"

print("Is file:", os.path.isfile(path))
print("Is dir :", os.path.isdir(path))

env = lmdb.open(path, readonly=True, lock=False, subdir=os.path.isdir(path))

with env.begin() as txn:
    for key, val in txn.cursor():
        print("KEY:", key)
        item = pickle.loads(val)
        print("FIELDS:", item.keys())
        break
import lmdb, pickle, numpy as np
with env.begin() as txn:
    rec = pickle.loads(txn.get(b"0"))
print("pocket_coord shape:", np.array(rec["pocket_coordinates"]).shape)
print("ligand shape:", np.array(rec["pocket_coord"]).shape if "pocket_coord" in rec else "no ligand")

PY
