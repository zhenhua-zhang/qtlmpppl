#!/bin/bash


export OUTDATED_IGNORE=1
module load Python/3.6.3-foss-2015b

source ../.env/bin/activate
pjdir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/wp_hiv_ccr5

# 一般而言大概需要25G内存

# hiv ccr5 expression percentage
python boxplot.py \
    -p ${pjdir}/outputs/age_gender/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_phtp.tsv \
    -c ${pjdir}/outputs/age_gender/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_cvrt.tsv \
    -g ${pjdir}/inputs/genotypes/200hiv/dosage/MatrixEQTL/200HIV_dosage.gz \
    -i ${pjdir}/inputs/genotypes/200hiv/dosage/MatrixEQTL/200HIV_variantInfo.gz \
    -P Pgated_CCR5P_mTreg_log10 \
    -G rs60939770 rs10897605 rs3176953 rs73934105 rs66942281 rs11574435 rs59440261 rs6441975 rs62246129 rs3087253 rs2373226 rs2213290 rs1001007 rs1015164 rs4317138 \
    -o ${pjdir}/outputs/age_gender/hiv.ccr5_exp.pc/qtlmapping/boxplot

# hiv ccr5 expression geometry mean
python boxplot.py \
    -p ${pjdir}/outputs/age_gender/hiv.ccr5_exp.gm/preprocess/hiv.ccr5_exp.gm.proc_phtp.tsv \
    -c ${pjdir}/outputs/age_gender/hiv.ccr5_exp.gm/preprocess/hiv.ccr5_exp.gm.proc_cvrt.tsv \
    -g ${pjdir}/inputs/genotypes/200hiv/dosage/MatrixEQTL/200HIV_dosage.gz \
    -i ${pjdir}/inputs/genotypes/200hiv/dosage/MatrixEQTL/200HIV_variantInfo.gz \
    -P GM_CCR5P_L GM_CCR5P_RApR7n_log10 GM_CCR5P_CD4 GM_CCR5P_CD8 GM_CCR5P_CM_CD8 GM_CCR5P_EM_CD8 GM_CCR5P_mTreg GM_CCR5P_RAnR7n GM_CCR5P_TEM "GM_CCR5P_RAnR7n_(EM)_CD8" "GM_CCR5P_RApR7n_(EM)_CD8" \
    -G rs10897605 rs3176953 rs73934105 rs66942281 rs11574435 rs60939770 rs59440261 rs6441975 rs62246129 rs3087253 rs2373226 rs2213290 rs1001007 rs1015164 rs4317138 \
    -o ${pjdir}/outputs/age_gender/hiv.ccr5_exp.gm/qtlmapping/boxplot


# bcg ccr5 expression cell propotion leading SNPs
# rs113010081,Pgated_CCR5P_CD8
# rs113010081,Pgated_CCR5P_L
# rs113341849,Pgated_CCR5P_CM_CD8_log10
# rs113341849,Pgated_CCR5P_EM_CD8
# rs113341849,Pgated_CCR5P_RAnR7n
# rs113341849,Pgated_CCR5P_RAnR7n_(EM)_CD8
# rs2628160,Pgated_CCR5P_M
# rs2798659,Pgated_CCR5P_Naive_CD8_log10
# rs35648233,Pgated_CCR5P_CD45
# rs75082326,Pgated_CCR5P_RAnR7p
# rs7826548,Pgated_CCR5P_nTreg_log10

# bcg ccr5 expression percentage
python boxplot.py \
    -p ${pjdir}/outputs/age_gender/bcg.ccr5_exp.pc/preprocess/bcg.ccr5_exp.pc.proc_phtp.tsv \
    -c ${pjdir}/outputs/age_gender/bcg.ccr5_exp.pc/preprocess/bcg.ccr5_exp.pc.proc_cvrt.tsv \
    -g ${pjdir}/inputs/genotypes/300bcg/dosage/MatrixEQTL/300BCG.dosage.tsv.gz \
    -i ${pjdir}/inputs/genotypes/300bcg/dosage/MatrixEQTL/300BCG.info.tsv.gz \
    -P "Pgated_CCR5P_CD8" "Pgated_CCR5P_L" "Pgated_CCR5P_CM_CD8_log10" "Pgated_CCR5P_EM_CD8" "Pgated_CCR5P_RAnR7n" "Pgated_CCR5P_RAnR7n_(EM)_CD8" "Pgated_CCR5P_M" "Pgated_CCR5P_Naive_CD8_log10" "Pgated_CCR5P_CD45" "Pgated_CCR5P_RAnR7p" "Pgated_CCR5P_nTreg_log10" \
    -G rs113010081 rs113341849 rs2628160 rs2798659 rs35648233 rs75082326 rs7826548  \
    -o ${pjdir}/outputs/age_gender/bcg.ccr5_exp.pc/qtlmapping/boxplot


# bcg ccr5 expression IFM leading SNPs
# rs113010081,GM_CCR5P_CD4
# rs113010081,GM_CCR5P_RAnR7n
# rs113010081,GM_CCR5P_RAnR7n_(EM)_CD8
# rs113010081,GM_CCR5P_TEM
# rs113341849,GM_CCR5P_CD8
# rs113341849,GM_CCR5P_mTreg
# rs113341849,GM_CCR5P_RApR7n
# rs113341849,GM_CCR5P_RApR7n_(EM)_CD8
# rs1441896,GM_CCR5P_CM_CD8
# rs17285324,GM_CCR5P_RAnR7p
# rs2117457,GM_CCR5P_M
# rs4775471,GM_CCR5P_RApR7p_log10
# rs9522137,GM_CCR5P_nTreg

# bcg ccr5 expression geometry mean
python boxplot.py \
    -p ${pjdir}/outputs/age_gender/bcg.ccr5_exp.gm/preprocess/bcg.ccr5_exp.gm.proc_phtp.tsv \
    -c ${pjdir}/outputs/age_gender/bcg.ccr5_exp.gm/preprocess/bcg.ccr5_exp.gm.proc_cvrt.tsv \
    -g ${pjdir}/inputs/genotypes/300bcg/dosage/MatrixEQTL/300BCG.dosage.tsv.gz \
    -i ${pjdir}/inputs/genotypes/300bcg/dosage/MatrixEQTL/300BCG.info.tsv.gz \
    -P "GM_CCR5P_CD4" "GM_CCR5P_RAnR7n" "GM_CCR5P_RAnR7n_(EM)_CD8" "GM_CCR5P_TEM" "GM_CCR5P_CD8" "GM_CCR5P_mTreg" "GM_CCR5P_RApR7n" "GM_CCR5P_RApR7n_(EM)_CD8" "GM_CCR5P_CM_CD8" "GM_CCR5P_RAnR7p" "GM_CCR5P_M" "GM_CCR5P_RApR7p_log10" "GM_CCR5P_nTreg" \
    -G rs113010081 rs113341849 rs1441896 rs17285324 rs2117457 rs4775471 rs9522137 \
    -o ${pjdir}/outputs/age_gender/bcg.ccr5_exp.gm/qtlmapping/boxplot

# python boxplot.py -p ../../outputs/age_gender/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_phtp.tsv -c ../../outputs/age_gender/hiv.ccr5_exp.pc/preprocess/hiv.ccr5_exp.pc.proc_cvrt.tsv -i ../../inputs/genotypes/200hiv/dosage/MatrixEQTL/per_chrom/chr3_variantInfo.gz -g ../../inputs/genotypes/200hiv/dosage/MatrixEQTL/per_chrom/chr3_dosage.gz -P Pgated_CCR5P_mTreg_log10 -G rs60939770 rs1015164
