#!/usr/bin/env Rscript
if(!require(dplyr)) install.packages("dplyr")
if(!require(stringi)) install.packages("stringi")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(reshape2)) install.packages("reshape2")
if(!require(data.table)) install.packages("data.table")

setwd("~/Documents/projects/wp_hiv_ccr5/outputs/age_gender_ivrk/integration/basicStatistics")
fp <- "~/Documents/projects/wp_hiv_ccr5/inputs/phenotypes/ccr5_expression.Dataset_3_Aug.newid.rmotl.csv"
df <- fread(fp, data.table=FALSE)

cnm <- names(df)
mtcols <- c("Cohort", "Sex", 'Age_at_visit')
mscols <- sapply(cnm, FUN=function(e) {
  (stri_startswith(e, fixed="Pgated") || stri_startswith(e, fixed="GM")) &&
    !stri_endswith(e, fixed="_G") && !stri_endswith(e, fixed="_RAp25p") || e %in% mtcols}, simplify=T)

df <- df[mscols]
df <- reshape2::melt(df, id=mtcols)

rnmap <- c(
  "CD45"= "Leukocytes",
  "M"="Monocytes",
  "L"="Lymphocytes",
  "CD4"="CD4+",
  "RApR7p"="CD4+Naive",
  "RAnR7p"="CD4+CM",
  "mTreg"="CD4+mTreg",
  "nTreg"="CD4+nTreg",
  "TEM"="CD4+TEM",
  "RAnR7n"="CD4+EM",
  "RApR7n"="CD4+TEMRA",
  "CD8"="CD8+",
  "Naive_CD8"="CD8+Naive",
  "CM_CD8"="CD8+CM",
  "EM_CD8"="CD8+TEM",
  "RAnR7n_EM_CD8"="CD8+EM",
  "RApR7n_EM_CD8"="CD8+TEMRA")

df <- df %>%
  mutate(Method=factor(ifelse(stri_startswith(variable, fixed="Pgated"), "CP", "MFI"), levels=c("MFI", "CP")),
         value=log2(value),
         variable=rnmap[stri_replace(variable, regex="(GM|Pgated)_CCR5P_", replacement = "")])

g <- ggplot(data=df) + theme_bw() +
  geom_boxplot(aes(y=reorder(variable, value, FUN=median), x=value, color=Cohort), outlier.size=0.5) +
  facet_wrap(Method~Cohort, scales="free") +
  labs(x="Log2 of measurements (MFI or CP)", y="Cell subpopulations")

ggsave("./MFIandCPPerCellSubpopulation.pdf", width=10, height=10)
ggsave("./MFIandCPPerCellSubpopulation.png", width=10, height=10)
