#!/bin/bash

host_name=$(hostname)

if [[ ${host_name} =~ calculon ]] || [[ ${host_name} =~ umcg-node ]]; then
    module purge
    module load Python/3.6.3-foss-2015b
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

#
## Convert excel sheet into tsv or csv
#
python ./src/convert_sheet_into_csv.py \
    --execle-file ${ipdir}/rawdata/Dataset_8_april.xlsx \
    --output-pre ${ipdir}/phenotypes/ \
    --sheet-index 'Dataset 5'

python ./src/convert_sheet_into_csv.py \
    --execle-file ${ipdir}/rawdata/Dataset_9_July.xlsx \
    --output-pre ${ipdir}/phenotypes/ccr5_exp \
    --sheet-index 'Dataset_9_July'

#
## You have to convert ID, as that the ID for each sample are different in genotypes and phenotypes
#

#
## Preprocess the phenotype and covariates
#
# mkdir -p ${outputdir} # You have to make directory based on you config file.
python ./src/preprocess.py \
    --conf-file ./configs/age_gender_cmv/hiv.ccr5_exp.pc.json


#
## QTL mapping by MatrixEQTL
#
module purge
module load R/3.5.1-foss-2015b-bare
module list

Rscript ./src/matrix_eqtl.r \
    --gntp-file ${ipdir}/genotype/dosage/MatrixEQTL/200HIV_dosage.gz \
    --gntp-info-file ${ipdir}/genotype/dosage/MatrixEQTL/200HIV_variantInfo.gz \
    --phtp-file ${opdir}/age_gender_cmv/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_phtp.tsv \
    --cvrt-file ${opdir}/age_gender_cmv/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_cvrt.tsv \
    --save-pref ${opdir}/age_gender_cmv/hiv.ccr5_exp.pc/qtlmapping/
