#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80
#SBATCH --time=24:00:00
#SBATCH -o logs/drugclip_test_DUDE_%j.out
#SBATCH -e logs/drugclip_test_DUDE_%j.err

module purge
module load cuda/11.8

source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate
cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome

python pocket_size.py

