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
python ./convert_sheet_into_csv.py \
    --execle-file ${ipdir}/rawdata/Dataset_8_april.xlsx \
    --output-pre ${ipdir}/phenotypes/ \
    --sheet-index 'Dataset 5'
