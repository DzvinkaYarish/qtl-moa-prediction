#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="borzoi_inference"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesla:1
#SBATCH -w falcon1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load cuda/11.7.0
module load cudnn/8.7.0.84-11.8

module load any/python/3.8.3-conda

conda activate borzoi_py39

${HOME}/.conda/envs/borzoi_py39/bin/python python_scripts/get_preds_borzoi.py \
${HOME}/qtl_labeling/moa_data/full_dataset_with_labeled_eqtls_and_negatives_for_borzoi.csv \
${HOME}/qtl_labeling/moa_data/borzoi_features_with_labeled_eqtls_and_negatives.csv \

