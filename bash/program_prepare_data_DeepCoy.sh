#!/bin/bash
#SBATCH --job-name=decoys
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o logs/preparedata_%j.out
#SBATCH -e logs/preparedata_%j.err

module purge
module load cuda/11.8
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

#version 1 fichier
#python deepcoy/data/prepare_data.py \
#  --data_path data/ValidationBiolip/1A4M/actives_clean.ism \
#  --dataset_name 1A4M \
#  --save_dir data/ValidationBiolip/1A4M/

#Version générale
for dir in data/ValidationEcoli_sansKiKd/*; do
    pdb=$(basename "$dir")
    
    input="$dir/actives_final.ism"
    
    # skip si fichier absent ou vide
    [ -s "$input" ] || { echo "⚠️ skip $pdb (pas de fichier)"; continue; }
    
    echo "🚀 Traitement $pdb"
    
    python deepcoy/data/prepare_data.py \
        --data_path "$input" \
        --dataset_name "$pdb" \
        --save_dir "$dir/"
done
