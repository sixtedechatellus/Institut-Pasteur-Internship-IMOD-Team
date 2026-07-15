#!/bin/bash
#SBATCH --job-name=encode_pockets_pcba
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH -o logs/encode_pockets_pcba_%j.out
#SBATCH -e logs/encode_pockets_pcba_%j.err

module purge
module load cuda/12.2
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate
cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

python unimol/data/encode_pocket.py \
  --user-dir unimol \
  --task drugclip \
  --arch drugclip \
  --results-path ./data/lit_pcba/lit_pcba \
  --batch-size 32 \
  --num-workers 8 \
  --fp16 \
  ./data/lit_pcba/lit_pcba
