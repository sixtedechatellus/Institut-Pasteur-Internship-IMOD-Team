#!/bin/bash

#SBATCH --job-name=benchmarkValidation
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:A40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/benchmark_%x_%A_%a.out
#SBATCH -e logs/benchmark_%x_%A_%a.err

# ========================
# Argument dataset
# ========================
if [ -z "$1" ]; then
  echo "Usage: sbatch $0 <nom_dataset>"
  exit 1
fi

DATASET=$1
FOLD=$SLURM_ARRAY_TASK_ID

BASE_DIR="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome"

module load cuda/11.8
source ${BASE_DIR}/bin/activate

cd ${BASE_DIR}
export PYTHONPATH=$PYTHONPATH:$(pwd)

results_path="${BASE_DIR}/résultats/${DATASET}"
batch_size=64
weight_path="./data/model_weights/6_folds/fold_${FOLD}.pt"
use_folds=True
TASK="BIOLIP"
data_path="${BASE_DIR}/data/${DATASET}/"

echo "[🚀] Début fold ${FOLD} - $(date)"

CUDA_VISIBLE_DEVICES="0" python ./unimol/test.py --user-dir ./unimol \
       $data_path \
       --results-path $results_path \
       --num-workers 8 --ddp-backend=c10d --batch-size $batch_size \
       --task drugclip --loss in_batch_softmax --arch drugclip \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256 --seed 1 \
       --use-folds $use_folds \
       --path $weight_path \
       --log-interval 100 --log-format simple \
       --max-pocket-atoms 511 \
       --test-task $TASK

echo "[✔] Fin fold ${FOLD} - $(date)"

# ========================
# Extraction automatique
# ========================

LOG_OUT="logs/benchmark_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${FOLD}.out"
CSV_OUT="résultats/${DATASET}_fold${FOLD}_6.csv"

echo "[📊] Extraction résultats fold ${FOLD}"

python extract_results2.py ${LOG_OUT} ${CSV_OUT}
