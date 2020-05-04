#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
#

## First: seperate dataframes for HIV and BCG cohors

HIV <- subset(Dataset_8_april, Cohort == "200HIV")
BCG <- subset(Dataset_8_april, Cohort == "300BCG")

## Second: select only the necessary variables

HIV2 <- HIV[, c(32:138, 5, 4)]
BCG2 <- BCG[, c(32:138, 5, 4)]

# Variable sex should be numeric

HIV2$Sex <- as.factor(HIV2$Sex)
levels(HIV2$Sex) <- c("0", "1")
HIV2$Sex <- as.numeric(HIV2$Sex)

BCG2$Sex <- as.factor(BCG2$Sex)
levels(BCG2$Sex) <- c("0", "1")
BCG2$Sex <- as.numeric(BCG2$Sex)

## Third: omit dataframes to remove NAs

HIV3 <- na.omit(HIV2)
BCG3 <- na.omit(BCG2)

## Fourth: script Hua for extracting correlation coefficients

library(MASS)
library(ppcor)

# merge pcor.test()$estimate results
merge_pctt <- function(data, ccr5_idx, bimk_idx, cvar_idx, mthd="spearman") {
    # Explanation of the arguments
    # data: the data.frame to be processed
    # ccr5_idx: the index of CCR5 expression measurements
    # bink_idx: the index of biomarks expressoin(?) measurements
    # cvar_idx: the index of covariates
    # mthd: the method to be used in pcor.test()

    ccr5_dtfm <- data[ccr5_idx]  # CCR5 expression data
    bimk_dtfm <- data[bimk_idx]  # Levels of biomarks
    cvar_dtfm <- data[cvar_idx]  # Co-variates

    ccr5_name_pool <- colnames(ccr5_dtfm)
    bimk_name_pool <- colnames(bimk_dtfm)

    cor_dtfm <- NULL
    for (ccr5_name in ccr5_name_pool) {
        ccr5_exp_vec <- ccr5_dtfm[ccr5_name]
        cor_vec <- c()
        for (bimk_name in bimk_name_pool) {
            bimk_exp_vec <- bimk_dtfm[bimk_name]
            pctt <- pcor.test(ccr5_exp_vec, bimk_exp_vec, cvar_dtfm, mthd)
            pctt_est <- pctt$estimate
            cor_vec <- c(cor_vec, pctt_est)
        }
        if (is.null(cor_dtfm))
            cor_dtfm <- cor_vec
        else
            cor_dtfm <- cbind(cor_dtfm, cor_vec)
    }

    colnames(cor_dtfm) <- ccr5_name_pool
    rownames(cor_dtfm) <- bimk_name_pool

    return(cor_dtfm)
    # If you want handle the results from one function to a variable, you have
    # use the return key words not print()
}


## Fifth: script Hua with small adjustments for extracting p-values

merge_pctt_pvalue <- function(data, ccr5_idx, bimk_idx, cvar_idx, mthd="spearman") {
    # Explanation of the arguments
    # data: the data.frame to be processed
    # ccr5_idx: the index of CCR5 expression measurements
    # bink_idx: the index of biomarks expressoin(?) measurements
    # cvar_idx: the index of covariates
    # mthd: the method to be used in pcor.test()
    
    ccr5_dtfm <- data[ccr5_idx]  # CCR5 expression data
    bimk_dtfm <- data[bimk_idx]  # Levels of biomarks
    cvar_dtfm <- data[cvar_idx]  # Co-variates
    
    ccr5_name_pool <- colnames(ccr5_dtfm)
    bimk_name_pool <- colnames(bimk_dtfm)
    
    cor_dtfm <- NULL
    for (ccr5_name in ccr5_name_pool) {
        ccr5_exp_vec <- ccr5_dtfm[ccr5_name]
        cor_vec <- c()
        for (bimk_name in bimk_name_pool) {
            bimk_exp_vec <- bimk_dtfm[bimk_name]
            pctt <- pcor.test(ccr5_exp_vec, bimk_exp_vec, cvar_dtfm, mthd)
            pctt_est_p <- pctt$p.value
            cor_vec <- c(cor_vec, pctt_est_p)
        }
        if (is.null(cor_dtfm))
            cor_dtfm <- cor_vec
        else
            cor_dtfm <- cbind(cor_dtfm, cor_vec)
    }
    
    colnames(cor_dtfm) <- ccr5_name_pool
    rownames(cor_dtfm) <- bimk_name_pool
    
    return(cor_dtfm)
    # If you want handle the results from one function to a variable, you have
    # use the return key words not print()
}


# Check if it works:

pctt <- merge_pctt(HIV3, 1:38, 39:107, 108:109)
print(pctt)

pctt2 <- merge_pctt_pvalue(HIV3, 1:38, 39:107, 108:109)
print(pctt2)

# Looking at variable Pgated_CCR5P_TEM and CSF-1 the estimate is 0.04813298 and the corresponding p-value is 0.5051115893

pcor.test(x = HIV3$Pgated_CCR5P_TEM, y = HIV3$`CSF-1`, z = c(HIV3$Age_at_visit, HIV3$Sex), method = "spearman")
# The results for this test are:     
#estimate   p.value statistic   n gp   Method
#1 0.04991094 0.3249293  0.985626 392  1 spearman

# Thus: there are some differences. I think it has something to do with the fact that we want to correct for two variables instead of one. Because if i only correct for age, it works perfectly:

pctt3 <- merge_pctt(HIV3, 1:38, 39:107, 108)
pctt4 <- merge_pctt_pvalue(HIV3, 1:38, 39:107, 108)

# Looking at variable Pgated_CCR5P_TEM and CSF-1 the estimate is 0.04762850 and the corresponding p-value is 0.5084872147

pcor.test(x = HIV3$Pgated_CCR5P_TEM, y = HIV3$`CSF-1`, z = HIV3$Age_at_visit, method = "spearman")
# The results for this test are:
# estimate   p.value statistic   n gp   Method
# 1 0.0476285 0.5084872  0.662428 196  1 spearman







