#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=64G
#SBATCH -o logs/drugclip_test_DUDE_%j.out
#SBATCH -e logs/drugclip_test_DUDE_%j.err

module load cuda/11.8
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

python - <<'PY'
import lmdb, pickle, os

# 🔧  Chemin ABSOLU correct :
path = "/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip/1A4M/pocket.lmdb"

import lmdb, pickle, numpy as np
env = lmdb.open(path, subdir=False, readonly=True)
with env.begin() as txn:
    cur = txn.cursor()
    keys = list(cur.iternext(keys=True, values=False))
    print("nb d'entrées:", len(keys))
    _, val = next(iter(cur))
    d = pickle.loads(val)
    print("atoms:", len(d["pocket_atoms"]))
    print("coordinates shape:", np.array(d["coordinates"]).shape)
