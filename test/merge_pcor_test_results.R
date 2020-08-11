#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
#

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

    ccr5_dtfm <- data[, ccr5_idx]  # CCR5 expression data
    bimk_dtfm <- data[, bimk_idx]  # Levels of biomarks
    cvar_dtfm <- data[, cvar_idx]  # Co-variates

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
    # If you want to hand the results from one function to a variable, you have
    # use the return key words not print()
}


# A data.frame for testing
dtfm <- data.frame(
   ccr5_a=1:10, ccr5_b=1:10, ccr5_c=1:10, bimk_a=1:10, bimk_b=1:10,
   bimk_c=1:10, age=21:30, gender=c(1, 1, 0, 1, 1, 0, 1, 0, 0, 1)
)

pctt <- merge_pctt(dtfm, 1:3, 4:6, 7:8)
print(pctt)
