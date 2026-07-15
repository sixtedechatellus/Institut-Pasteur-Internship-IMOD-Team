#!/bin/bash
#SBATCH --job-name=encode_mols_pcba
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH -o logs/encode_mols_pcba_%j.out
#SBATCH -e logs/encode_mols_pcba_%j.err

module purge
module load cuda/12.2
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate
cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

# Lancer l'encodage
python unimol/encode_mols.py \
   --user-dir unimol \
   --task drugclip \
   --arch drugclip \
   --results-path ./data/lit_pcba/lit_pcba \
   --batch-size 64 \
   --num-workers 8 \
   --fp16 \
   --mol-path ./data/lit_pcba/lit_pcba \
   ./data/lit_pcba/lit_pcba
