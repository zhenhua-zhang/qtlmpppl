#!/bin/bash
sbop=age_gender_cmv
sbop=age_gender

mm=hiv.ccr5_exp.gm
mm=hiv.ccr5_exp.pc
mm=bcg.ccr5_exp.pc
#mm=bcg.ccr5_exp.gm

qos=regular

sbatch \
    --time 4:59:0 \
    --mem 80G \
    --cpus-per-task 1 \
    --job-name ${mm} \
    --qos ${qos} \
    --output %j-%u-${mm} <<EOF
#!/bin/bash

source /apps/modules/modules.bashrc

pjdir=\${HOME}/Documents/projects/wp_hiv_ccr5
ipdir=\${pjdir}/inputs
opdir=\${pjdir}/outputs
scdir=\${pjdir}/scripts
sbopdir=\${opdir}/${sbop}

host_name=\$(hostname)
if [[ \${host_name} =~ calculon ]] || [[ \${host_name} =~ umcg-node ]]; then
    module purge
    module load Python/3.6.3-foss-2015b
    source \${scdir}/.env/bin/activate
    module list
fi

set -eu -o pipefail

mkdir -p \${sbopdir}/${mm}/{preprocess,qtlmapping/{per_trait,manhattan,boxplot,qq,locuszoom},function_annotation/FUMA}

python \${scdir}/src/preprocess.py \
    --conf-file \${scdir}/configs/${sbop}/${mm}.json

if [[ \${host_name} =~ calculon ]] || [[ \${host_name} =~ umcg-node ]]; then
    module purge
    module load R/3.5.1-foss-2015b-bare
    module list
fi


Rscript \${scdir}/src/matrix_eqtl.r \
    --gntp-file \${ipdir}/genotypes/300bcg/dosage/MatrixEQTL/300BCG.dosage.tsv.gz \
    --gntp-info-file \${ipdir}/genotypes/300bcg/dosage/MatrixEQTL/300BCG.info.tsv.gz \
    --phtp-file \${opdir}/${sbop}/${mm}/preprocess/${mm}.proc_phtp.tsv \
    --cvrt-file \${opdir}/${sbop}/${mm}/preprocess/${mm}.proc_cvrt.tsv \
    --save-pref \${opdir}/${sbop}/${mm}/qtlmapping/
EOF
