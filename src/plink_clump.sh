#!/bin/bash
set -Eeu -o pipefail

# There are duplicated SNP ids in Chr2. One has to removed them manually.
# wkdir=~/Documents/projects/wp_hiv_ccr5/inputs/references/plink/1kg_phase3
# for chrom in chr{2,3}; do
#     if [[ ! -e ${wkdir}/EUR.${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.af0.01.vcf.gz ]]; then
#         bcftools view \
#             -i 'TYPE="snp" && INFO/EUR_AF>0.01' \
#             -S ${wkdir}/1kg_phase3-eur_id.txt \
#             -o ${wkdir}/EUR.${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.af0.01.vcf.gz \
#             -O z \
#             --threads 4 \
#             ${wkdir}/ALL.${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#     fi
# 
#     if [[ ! -e ${wkdir}/eur-1kg_p3-${chrom}.log ]]; then
#         plink --vcf ${wkdir}/EUR.${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.af0.01.vcf.gz \
#             --make-bed \
#             --out ${wkdir}/eur-1kg_p3-${chrom}
#     fi
# 
#     if [[ ! -e ${wkdir}/eur-1kg_p3-${chrom}-dup.log ]]; then
#         plink --bfile ${wkdir}/eur-1kg_p3-${chrom} \
#             --list-duplicate-vars ids-only suppress-first \
#             --out ${wkdir}/eur-1kg_p3-${chrom}-dup
#     fi
# 
#     if [[ $(head ${wkdir}/eur-1kg_p3-${chrom}-dup.dupvar) -ge 1 ]]; then
#         plink --bfile ${wkdir}/eur-1kg_p3-${chrom} \
#             --exclude suppress-first ${wkdir}/eur-1kg_p3-${chrom}-dup.dupvar \
#             --make-bed \
#             --out ${wkdir}/eur-1kg_p3-${chrom}-nodup
# 
#         for fmt in bim fam bed; do 
#             mv ${wkdir}/eur-1kg_p3-${chrom}-nodup.${fmt} ${wkdir}/eur-1kg_p3-${chrom}.${fmt} -f
#         done
#     else
#         rm ${wkdir}/eur-1kg_p3-${chrom}-dup* -f
#     fi
# done

pjdir=~/Documents/projects/wp_hiv_ccr5
clump_dir=${pjdir}/outputs/age_gender/integration/plink_clump/"$1"

do_clump=0
declare -A prefix_map=([gm]=GM_CCR5P [pc]=Pgated_CCR5P)
chrom=3
if [[ ${do_clump} -eq 1 ]]; then
    for cohort in bcg hiv; do
        for mthd in pc gm; do
            prodir=${pjdir}/outputs/age_gender/${cohort}.ccr5_exp.${mthd}
            for celltype in L M TEM CD{4,8,45} {CM,EM,Naive}_CD8 RA{n,p}R7n_EM_CD8 {m,n}Treg RA{n,p}R7{n,p}; do
                ssfile=${prodir}/qtlmapping/per_trait/${prefix_map[$mthd]}_${celltype}.csv.gz
                if [[ ! -e ${ssfile} ]]; then
                    ssfile=${ssfile/.csv.gz/_log10.csv.gz}
                    if [[ ! -e ${ssfile} ]]; then
                        ssfile=${ssfile/_log10/_ivrk}
                    fi
                fi

                sssfile=${clump_dir}/${cohort}-${mthd}-${prefix_map[$mthd]}_${celltype}-chr${chrom}.tsv
                if [[ ! -e ${sssfile} ]]; then
                    zcat "${ssfile}" | awk -F',' '{if(($4<5e-8 && $7=='$chrom')||$4=="pvalue") {print $1"\t"$4"\t"$6}}' > "${sssfile}"
                fi

                if [[ $(grep rs "${sssfile}" -c) -gt 1 ]]; then
                    clump_file=${sssfile/.tsv/-plink_clump}
                    plink --bfile ~/Documents/projects/wd_plink_data/1kg_phase3/eur-1kg_p3-chr$chrom \
                        --clump "${sssfile}" \
                        --clump-field pvalue \
                        --clump-snp-field snps \
                        --clump-p1 7.35e-10 \
                        --clump-p2 5e-8 \
                        --clump-r2 0.6 \
                        --clump-kb 500 \
                        --out "${clump_file}"
                else
                    rm -f "${sssfile}"
                fi

            done
        done
    done
fi

# Merge independent variants
miv=0
if [[ ${miv} -eq 1 ]]; then
    cd ${clump_dir} || exit
    (echo -e 'Cohort\tMeasurement\tCellType\tChrom\tIndexSNP\tPos\tPvalue\tClumpedSNPs'
    grep rs ./*.clumped \
        | sed -e 's/ \+/\t/g; s/://g; s/\(Pgated\|GM\)_CCR5P_\(.*\)-chr3-plink_clump.clumped/\2/g' \
        -e 's/hiv/HIV patients/g; s/bcg/Healthy controls/g' \
        -e 's/-pc-/\tCPCp\t/g; s/-gm-/\tMFI\t/g' \
        -e 's#./##g' \
        | cut -f1-4,6-8,15) \
        | awk -F '\t' '{print $5"\t"$4"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"\t"$8}' \
        | sort -k1,1g \
        > ./independent_variants.tsv
    cd ..
fi

# Prepare data for bubble plot
pp_bbp=0
if [[ ${pp_bbp} -eq 1 ]]; then
    (echo "#cohort,mthd,celltype,pre_clump,post_clump"
    for x in "${clump_dir}"/*.clumped; do
        pre_clump_frp=$(readlink -f "${x/-plink_clump.clumped/.tsv}")
        pst_clump_frp=$(readlink -f "${x}")
        meta_info=$(cut -f1-3 -d'-' --output-delimiter=',' <<<$(basename "${x}"))
        echo "$meta_info,${pre_clump_frp},${pst_clump_frp}"
    done ) > ${clump_dir}/bubble_plot-meta_file.csv

    python ${pjdir}/scripts/src/plink_clump-bubble_plot.py \
        -i ${clump_dir}/bubble_plot-meta_file.csv \
        -F png pdf \
        -o ${clump_dir}/clumped_snps
fi


# 需要先手动把所有clumped的文件组合到一起。
# while read -r line; do
#     cohort=$(cut -f1 -d, <<<"$line")
#     mthd=$(cut -f2 -d, <<<"$line")
#     celltype=$(cut -f3 -d, <<<"$line")
#     snp=$(cut -f6 -d, <<<"$line")
#     file="$cohort-$mthd*_CCR5P_$celltype.tsv"
#     (grep -w "$snp" "${file}" | cut -f1,3 -d$'\t' | tr '\t\n' ','; echo "$line")
# done < <(grep -v Cohort ../ccr5_exp_qtl-chr3-plink_clumped.tsv) \
#     > ../ccr5_exp_qtl-chr3-plink_clumped.csv
