#!/usr/bin/env Rscript

library(dplyr)
library(ggpubr)
library(ggplot2)
library(data.table)

# plot violin plot
violin_plot <- function(dttb, var_name, p_title) {
  g <- ggplot(data=dttb) + theme_bw() +
    geom_boxplot(aes_string(y=var_name, x="classes"), alpha=0.5) +
    geom_jitter(aes_string(y=var_name, x="classes"), alpha=0.5, size=1) +
    ylab(p_title) +
    scale_x_discrete(element_blank())
  return(g)
}

pjdir <- "~/Documents/projects/wp_hiv_reservoir"
setwd(paste0(pjdir, "/outputs/hiv_res/preprocess"))

cvrt_flnm <- paste0(pjdir, "/outputs/hiv_res/preprocess/hiv-reservoir.proc_cvrt.tsv")
cvrt_dttb <- fread(cvrt_flnm)

phtp_flnm <- paste0(pjdir, "/outputs/hiv_res/preprocess/hiv-reservoir.proc_phtp.tsv")
phtp_dttb <- fread(phtp_flnm)

name_vec <- c("V1"="gender", "V2"="age", "V3"="CD4_NADIR", "V4"="HIV_DURATION", "V5"="RNAHIV_CD4LOG", "V6"="DNAHIV_CD4LOG", "V7"="RNAvsDNA_CD4LOG")
dttb <- rbindlist(list(cvrt_dttb, phtp_dttb)) %>%
  select(starts_with("X")) %>%
  transpose() %>%
  rename_with(function(e) {return(name_vec[e])}) %>%
  mutate(classes=1)

title_vec <- c("age"="Age (years)",
               "HIV_DURATION"="HIV duration (years)",
               "CD4_NADIR"="CD4 nadir",
               "RNAHIV_CD4LOG"="CA HIV RNA",
               "DNAHIV_CD4LOG"="Total HIV DNA",
               "RNAvsDNA_CD4LOG"="RRD")

var_name <- "age"
g_age <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

var_name <- "HIV_DURATION"
g_hiv_du <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

var_name <- "CD4_NADIR"
g_cd4_na <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

var_name <- "RNAHIV_CD4LOG"
g_hiv_rna <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

var_name <- "DNAHIV_CD4LOG"
g_hiv_dna <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

var_name <- "RNAvsDNA_CD4LOG"
g_rna_vs_dna <- violin_plot(dttb, var_name, p_title=title_vec[var_name])

plotslist <- ggarrange(g_hiv_dna, g_hiv_rna, g_rna_vs_dna, g_age, g_hiv_du, g_cd4_na, ncol=6, nrow=1)
ggexport(plotslist=plotslist, filename="host_traits_boxplot.pdf", width = 560)
ggexport(plotslist=plotslist, filename="host_traits_boxplot.png", width = 560)
