#!/bin/bash
#!/bin/bash
#SBATCH --job-name=create_mols
#SBATCH --output=logs/create_mols_%j.out
#SBATCH --error=logs/create_mols_%j.err
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Dataset path
DATASET_PATH="data/Alphafold"

python script_create_mols_alphafold.py "$DATASET_PATH"
