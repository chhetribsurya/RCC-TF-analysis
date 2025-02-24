#!/usr/bin/bash

bedfile_list="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project/HEK293_bed"
script_path="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project/scripts"
workdir="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project"
#outdir="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project/HEK293_output"

#gwas_list="/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project/GWAS_RCC"
#ld_blocklist="2500000 2000000 1000000 500000"

#outputdir="${workdir}/sbatch_outdir" 
logdir="$workdir/logfiles"
mkdir -p $logdir

#source activate r-env
eval "$(conda shell.bash hook)"
conda activate r-env

#for bedfile in ${bedfile_list}/*.bed; do
for bedfile in ${bedfile_list}/redo_bed/*.bed; do
    #for gwas in ${gwas_list}/*.tsv; do
        #for ldblock in ${ld_blocklist};do
            bedname=$(basename $bedfile)
            sbatch -p defq -n 15 --mem 60G -t 10:00:00 -o $logdir/${bedname}.out -e $logdir/${bedname}.err -J "${bedname}-RCCanalysis" --wrap="Rscript --vanilla ${script_path}/rcc_script.R ${bedfile} ${workdir}";
             #Rscript ${script_path}/rcc_script.R ${bedfile} ${workdir};
        #done
    #done
done
