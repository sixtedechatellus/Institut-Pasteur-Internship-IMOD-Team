#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o logs/decoys_%j.out
#SBATCH -e logs/decoys_%j.err

# Program to clean the actives SMILES from PDB protein/actives complexes for later DeepCoy use. See clean_smiles_files4.py to specify the folder
module purge
module load cuda/11.8
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

python clean_smiles_files4.py
