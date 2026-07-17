#!/bin/bash
#SBATCH --job-name=drugclip_test_DUDE
#SBATCH --partition=dedicatedgpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o logs/genfiles_%j.out
#SBATCH -e logs/genfiles_%j.err

module purge
module load cuda/11.8
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

python programme_extract_file_multiple_ligands.py


