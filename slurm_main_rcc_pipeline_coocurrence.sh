#!/bin/bash

#SBATCH --job-name=TF-Coocur-Analysis
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=1-936%110
#SBATCH --time=3:00:00
#SBATCH --mem=110G
#SBATCH --cpus-per-task=35
#SBATCH --partition=defq

#previous-array=1-1515%200

# Load environment
eval "$(conda shell.bash hook)"
conda activate /home/schhetr1/anaconda3/envs/r-env
#conda activate r-env

# Define the working directory
workdir="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project"

# Define the output directory
outdir="/data/abattle4/surya/datasets/for_diptavo/RCC_TF_Project/HEK293_expanded_bed_updatedThresh1e4"

mkdir -p $outdir

# Directory where bedfiles are stored
#bedfile_dir="${workdir}/HEK293_expanded_bed"
bedfile_dir="${workdir}/HEK293_expanded_bed_updatedThresh1e4"

# CSV file with bedfile pairs
#csv_file="${workdir}/test_HEK293_BATCH_output/AllminusEncode_1515_TF_cooccurrence_list.csv"
csv_file="${workdir}/test_HEK293_BATCH_output/UpdatedthreshminusPreviousthresh_936_TF_cooccurrence_list.csv"
#csv_file="${workdir}/test_HEK293_BATCH_output/AllminusEncode_1515_TF_cooccurrence_list_test.csv"

# Extract the bedfile names for the current array task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $csv_file)
bedfile_1=$(echo $line | cut -d ',' -f 1)
bedfile_2=$(echo $line | cut -d ',' -f 2)

# Full paths to bedfiles
bedfile_1_path="${bedfile_dir}/${bedfile_1}"
bedfile_2_path="${bedfile_dir}/${bedfile_2}"

# Print the bedfile paths to the Slurm log file
echo "Processing bedfile 1: $bedfile_1_path"
echo "Processing bedfile 2: $bedfile_2_path"

# Script path
script_path="./rcc_script_coccur_betaSE.R"

# Execute the R script with the bedfile paths and workdir
Rscript --vanilla "$script_path" "$bedfile_1_path" "$bedfile_2_path" "$workdir" "$outdir"
