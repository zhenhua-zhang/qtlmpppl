#!/bin/bash

host_name=$(hostname)

if [[ ${host_name} =~ calculon ]] || [[ ${host_name} =~ umcg-node ]]; then
    module purge
    module load Python
    module list
fi

if [ -d .env ]; then
    source ./.env/bin/activate
fi

pjdir=~/Documents/projects/wp_hiv_ccr5
ipdir=${pjdir}/inputs
opdir=${pjdir}/outputs
scdir=${pjdir}/scripts
tmdir=${pjdir}/temps

# Convert excel sheet into tsv or csv
python ./src/convert_sheet_into_csv.py \
    --execle-file ${ipdir}/rawdata/Dataset_8_april.xlsx \
    --output-pre ${ipdir}/phenotypes/ \
    --sheet-index 'Dataset 5'

# Preprocess the phenotype and covariates
python ./src/preprocess.py \
    --conf-file ./configs/hiv.ccr5_exp.pc.json


# QTL mapping by MatrixEQTL
module purge
module load R/3.5.1
module list

Rscript ./src/matrix_eqtl.r \
    --gntp-file ${ipdir}/genotype/dosage/MatrixEQTL/200HIV_dosage.gz \
    --gntp-info-file ${ipdir}/genotype/dosage/MatrixEQTL/200HIV_variantInfo.gz \
    --phtp-file ${opdir}/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_phtp.tsv \
    --cvrt-file ${opdir}/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_cvrt.tsv \
    --save-pref ${opdir}/hiv.ccr5_exp.pc/qtlmapping/
