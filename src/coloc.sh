#!/bin/bash
pjdir=~/Documents/projects/wp_hiv_ccr5
mafp=${pjdir}/outputs/age_gender/integration/maf/200hiv_snps_maf.tsv

vlss_dir=~/Documents/projects/wd_gwas_summary_statistics/PJ_McLaren_etal-polymorphisms-PNAS-2015/

declare -A prefix_map=([gm]=GM_CCR5P [pc]=Pgated_CCR5P)
declare -A mthd_map=([gm]="Geometry mean" [pc]="Cell proportion")
declare -A ssize_map=([hiv]=212 [bcg]=304)

rm -f ${pjdir}/outputs/age_gender/integration/coloc/{ccr5-vs-viral_load-per_cell_type.txt,comparison_cohort-per_cell_type.txt,comparison_measurements-per_cell_type.txt}

header=CellType,Measurement,Chrom,WithSigHits,PP0,PP1,PP2,PP3,PP4
if [[ ! -e ${pjdir}/outputs/age_gender/integration/coloc/ccr5-vs-viral_load-per_cell_type.txt ]]; then
    echo -e "${header}" > ${pjdir}/outputs/age_gender/integration/coloc/ccr5-vs-viral_load-per_cell_type.txt
fi

if [[ ! -e ${pjdir}/outputs/age_gender/integration/coloc/comparison_cohort-per_cell_type.txt ]]; then
    echo -e "${header}" > ${pjdir}/outputs/age_gender/integration/coloc/comparison_cohort-per_cell_type.txt
fi

if [[ ! -e ${pjdir}/outputs/age_gender/integration/coloc/comparison_measurements-per_cell_type.txt ]]; then
    echo -e "${header}" > ${pjdir}/outputs/age_gender/integration/coloc/comparison_measurements-per_cell_type.txt
fi

for chrom in 2 3 4 13; do
    for cell_type in L M TEM CD{4,8,45} {CM,EM,Naive}_CD8 RA{n,p}R7n_EM_CD8 {m,n}Treg RA{n,p}R7{n,p}; do
        vlss_ss=${vlss_dir}/meta_eur.${chrom}.meta
        if [[ (($chrom -eq 2 || $chrom -eq 4)) && $cell_type != "L" ]]; then continue; fi
        if [[ $chrom -eq 13 && $cell_type != "RAnR7p" ]]; then continue; fi

        for ccr5_dir in {bcg,hiv}.ccr5_exp.{pc,gm}; do
            if [[ ${chrom} -ne 3 ]]; then continue; fi
            pplt=${ccr5_dir%%.*}
            mthd=${ccr5_dir##*.}
            prefix=${prefix_map[${mthd}]}
            ccr5_ss=${pjdir}/outputs/age_gender/${ccr5_dir}/qtlmapping/per_trait/${prefix}_${cell_type}.csv.gz
            if [[ ! -e ${ccr5_ss} ]]; then
                ccr5_ss=${ccr5_ss/.csv.gz/_log10.csv.gz}
                if [[ ! -e ${ccr5_ss} ]]; then ccr5_ss=${ccr5_ss/_log10.csv.gz/_ivrk.csv.gz}; fi
            fi
            ccr5_ok=$(zcat "${ccr5_ss}"| awk -F',' '{if($4<5e-8 && $7=='$chrom'){print}}'| grep . -cm 1)
            vlss_ok=$(sed 's/ \+/,/g; s/^,//g' "${vlss_ss}" | awk -F',' '{if($7<5e-8 && $1=='$chrom'){print}}'| grep . -cm 1)

            [[ ${ccr5_ok} -ge 1 || ${vlss_ok} -ge 1 ]] \
                && echo -e "${cell_type},${pplt}-${mthd_map[${mthd}]},${chrom},ccr5(${ccr5_ok})|viral(${vlss_ok})",$(${pjdir}/scripts/src/coloc.r -a "${ccr5_ss}" -b "${vlss_ss}" -B SNP,P,CHR,BP -m ${mafp} -z ${ssize_map[${pplt}]},6057 -c ${chrom} | sed '/PP/d; s/ \+/,/g; s/^,//g; s/,$//g')
        done >> ${pjdir}/outputs/age_gender/integration/coloc/ccr5-vs-viral_load-per_cell_type.txt

        for mthd in pc gm; do
            case_dir=${pjdir}/outputs/age_gender/hiv.ccr5_exp.${mthd}/qtlmapping/per_trait
            ctrl_dir=${pjdir}/outputs/age_gender/bcg.ccr5_exp.${mthd}/qtlmapping/per_trait

            prefix=${prefix_map[${mthd}]}
            case_ss="${case_dir}/${prefix}_${cell_type}.csv.gz"
            if [[ ! -e ${case_ss} ]]; then
                case_ss=${case_ss/.csv.gz/_log10.csv.gz}
                if [[ ! -e ${case_ss} ]]; then case_ss=${case_ss/_log10.csv.gz/_ivrk.csv.gz}; fi
            fi
            case_ok=$(zcat "${case_ss}"| awk -F',' '{if($4<5e-8 && $7=='$chrom'){print}}'| grep . -cm 1)

            ctrl_ss="${ctrl_dir}/${prefix}_${cell_type}.csv.gz"
            if [[ ! -e ${ctrl_ss} ]]; then
                ctrl_ss=${ctrl_ss/.csv.gz/_log10.csv.gz}
                if [[ ! -e ${ctrl_ss} ]]; then ctrl_ss=${ctrl_ss/_log10.csv.gz/_ivrk.csv.gz}; fi
            fi
            ctrl_ok=$(zcat "${ctrl_ss}"| awk -F',' '{if($4<5e-8 && $7=='$chrom'){print}}'| grep . -cm 1)

            [[ ${case_ok} -ge 1 || ${ctrl_ok} -ge 1 ]] \
                && echo -e "${cell_type},${mthd_map[${mthd}]},${chrom},hiv(${case_ok})|bcg(${ctrl_ok})",$(${pjdir}/scripts/src/coloc.r -a "${ctrl_ss}" -b "${case_ss}" -m ${mafp} -z 304,212 -c ${chrom} | sed '/PP/d; s/ \+/,/g; s/^,//g; s/,$//g')
        done >> ${pjdir}/outputs/age_gender/integration/coloc/comparison_cohort-per_cell_type.txt

        for cohort in bcg hiv; do
            pc_dir=${pjdir}/outputs/age_gender/${cohort}.ccr5_exp.pc/qtlmapping/per_trait
            gm_dir=${pjdir}/outputs/age_gender/${cohort}.ccr5_exp.gm/qtlmapping/per_trait

            pc_ss="${pc_dir}/Pgated_CCR5P_${cell_type}.csv.gz"
            if [[ ! -e ${pc_ss} ]]; then
                pc_ss=${pc_ss/.csv.gz/_log10.csv.gz}
                if [[ ! -e ${pc_ss} ]]; then pc_ss=${pc_ss/_log10.csv.gz/_ivrk.csv.gz}; fi
            fi
            pc_ok=$(zcat "${pc_ss}"| awk -F',' '{if($4<5e-8 && $7=='$chrom'){print}}'| grep . -cm 1)

            gm_ss="${gm_dir}/GM_CCR5P_${cell_type}.csv.gz"
            if [[ ! -e ${gm_ss} ]]; then
                gm_ss=${gm_ss/.csv.gz/_log10.csv.gz}
                if [[ ! -e ${gm_ss} ]]; then gm_ss=${gm_ss/_log10.csv.gz/_ivrk.csv.gz}; fi
            fi
            gm_ok=$(zcat "${gm_ss}"| awk -F',' '{if($4<5e-8 && $7=='$chrom'){print}}'| grep . -cm 1)

            [[ ${pc_ok} -ge 1 || ${gm_ok} -ge 1 ]] \
                && echo -e "${cell_type},${cohort},${chrom},cp(${pc_ok})|gm(${gm_ok})",$(${pjdir}/scripts/src/coloc.r -a "${pc_ss}" -b "${gm_ss}" -m ${mafp} -z 304,212 -c ${chrom} | sed '/PP/d; s/ \+/,/g; s/^,//g; s/,$//g')
        done >> ${pjdir}/outputs/age_gender/integration/coloc/comparison_measurements-per_cell_type.txt
    done
done

unset mthd_map prefix_map ssize_map

Rscript ${pjdir}/scripts/src/coloc-buble_plot.r "$(pwd)"
