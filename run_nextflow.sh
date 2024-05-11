#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="moa"
#SBATCH --partition=main
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load any/jdk/1.8.0_265
module load any/singularity/3.7.3
module load squashfs/4.4
module load nextflow

export root=/gpfs/space/home/dzvenymy

#-stub-run
nextflow main.nf  --in_file "${root}/qtl_labeling/moa_data/full_dataset_with_labeled_eqtls_and_negatives.csv"  --HOME $root --cell_type 156 \
                --out_dir "${root}/Thesis/nextflow_output_full_dataset_cropped_enformer"

#CSV_DIR="/gpfs/space/home/dzvenymy/qtl_labeling/cQTL_data"
#
#for csv_file in "$CSV_DIR"/*.csv; do
#
#  if [ -e "${root}/qtl_labeling/cQTL_data/Enformer_crop/${csv_file%.*}_enformer_preds.csv" ]; then
#    echo "File exists"
#  else
#    # Check if the file is a regular file (and not an empty pattern)
#    if [ -f "$csv_file" ]; then
#        echo "Processing $csv_file"
#        # Add your processing commands here
#        # For example, you might want to count lines in each CSV
#        nextflow enformer_main.nf  --HOME $root --cell_type 156 \
#                --out_dir "${root}/qtl_labeling/cQTL_data/Enformer_crop/" --in_file $csv_file
#    fi
#  fi
#done


