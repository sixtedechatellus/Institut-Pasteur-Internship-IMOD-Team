#!/bin/bash
#SBATCH --job-name=benchmarkValidationEcoli
#SBATCH --partition=gpu
#SBATCH --qos=fast
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=64G
#SBATCH -o logs/benchmark_ValidationEcoli_fold5_%j.out
#SBATCH -e logs/benchmark_ValidationEcoli_fold5_%j.err

module load cuda/11.8
source /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/bin/activate

cd /pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome
export PYTHONPATH=$PYTHONPATH:$(pwd)

results_path="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/résultats/ValidationEcoli"
batch_size=64
weight_path="./data/model_weights/6_folds/fold_5.pt"
use_folds=True
TASK="BIOLIP"
data_path="/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationEcoli/"  # à vérifier

echo "[🚀] Début du test $(date)"

CUDA_VISIBLE_DEVICES="0" python ./unimol/test.py --user-dir ./unimol \
	 $data_path  \
       --results-path $results_path \
       --num-workers 8 --ddp-backend=c10d --batch-size $batch_size \
       --task drugclip --loss in_batch_softmax --arch drugclip  \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256  --seed 1 \
       --use-folds $use_folds \
       --path $weight_path \
       --log-interval 100 --log-format simple \
       --max-pocket-atoms 511 \
       --test-task $TASK \

echo "[✔] Fin du test $(date)"

# lignes dans test_dude pour lancer
#data_path = getattr(self.args, "data", None)
#        if data_path is None or not os.path.exists(data_path):
#            data_path = "./data/DUD-E"

#        print(f"[INFO] >>> Using dataset at: {data_path}")
#        targets = os.listdir(data_path) 

