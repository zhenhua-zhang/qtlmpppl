#!/usr/bin/env Rscript
library(data.table)

setwd('~/Documents/projects/wp_hiv_reservoir')

phtp <- fread('./inputs/phenotype/hiv_reservoir/20190524_HIVreservoir_GENT_withRNADNARatio.tsv', header = TRUE, data.table = FALSE)
row.names(phtp) <- phtp$id
phtp <- phtp[, !colnames(phtp) %in% c('id', 'COMMENTS_HIVDNA')]

cvrt <- fread('./inputs/phenotype/trait/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_cARTDuration_CD4Latest_20200212.csv', header = TRUE, data.table = FALSE)
row.names(cvrt) <- cvrt$id
cvrt <- cvrt[, !colnames(cvrt) %in% c('id')]

gntp <- fread('./misc/test/test_snps-dosage.txt', header = TRUE, data.table = FALSE)
row.names(gntp) <- gntp$id
gntp <- gntp[, !colnames(gntp) %in% c('id')]
gntp <- as.data.frame(t(gntp))

can_ind <- intersect(intersect(rownames(phtp), rownames(cvrt)), rownames(gntp))
lmdf <- cbind(gntp[can_ind, ], cvrt[can_ind, ],  phtp[can_ind, ])

m <- lm('DNAHIV_CD4LOG ~ RNAHIV_CD4LOG + rs7113204 + age + gender + CD4_NADIR + HIV_DURATION', data = lmdf)
summary(m)
