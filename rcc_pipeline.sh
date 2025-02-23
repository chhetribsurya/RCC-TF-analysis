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

#https://unix.stackexchange.com/questions/490649/pairwise-combinations-of-filenames
#https://unix.stackexchange.com/questions/338635/create-combinations-of-elements-from-a-vector-to-give-as-input-in-a-program
#declare -a encode_ids=(ENCFF002CDP ENCFF002COQ ENCFF002DAJ ENCFF002DCM)
#[schhetr1@login03 RCC_TF_Project]$ for i in ${encode_ids[@]}; do for j in ${encode_ids[@]}; do if [ "$i" \< "$j" ]; then echo "Pairs $i and $j";fi;done;done
#Pairs ENCFF002CDP and ENCFF002COQ
#Pairs ENCFF002CDP and ENCFF002DAJ
#Pairs ENCFF002CDP and ENCFF002DCM
#Pairs ENCFF002COQ and ENCFF002DAJ
#Pairs ENCFF002COQ and ENCFF002DCM
#Pairs ENCFF002DAJ and ENCFF002DCM
