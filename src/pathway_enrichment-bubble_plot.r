#!/usr/bin/env Rscript
library(dplyr)
library(stringi)
library(ggplot2)
library(data.table)

args = commandArgs(TRUE)

setwd(args[1])

lddt <- function(ptrn) {
    #' Load dataset.
    #'
    #' @param ptrn character pattern of files
    #'
    #' @return data.table
    
    files <- list.files(pattern = ptrn)
    
    whldt <- NULL
    for (fn in files) {
        bsfn <- stri_replace(fn, "", fixed = "-gene_sets.txt")
        fn_vec <- stri_split(bsfn, fixed = ".")[[1]]
        cht <- fn_vec[1]
        mmt <- fn_vec[3]
        ctp <- fn_vec[4]
        
        tmpdt <- fread(fn, verbose = FALSE)
        tmpdt[, `:=`(cohort = cht, measurement = mmt, celltype = ctp)]
        
        if (is.null(whldt)) { whldt <- tmpdt }
        else { whldt <- rbind(whldt, tmpdt) }
    }
    
    return(whldt)
}

main <- function() {
    # Load data and select category
    if (length(args) == 2) {
      cate <- args[2]
    } else {
      cate <- "GO_bp"
    }

    whldt <- lddt("*.txt") %>%
      filter(Category == cate, -log10(adjP) > 2)

    # Reorder the GeneSet by adjP
    whldt[order(GeneSet, -log10(adjP))]
    gene_sets <- unique(whldt$GeneSet)
    whldt$GeneSet <- factor(whldt$GeneSet, levels = gene_sets)
    y_ticks_labels <- stri_trans_totitle(stri_replace(gene_sets, "", fixed="GO_"))

    x_celltype_levels <- c("CD45", "L", "CD4", "mTreg", "TEM", "RApR7n", "RAnR7n", "CD8", "CM_CD8", "EM_CD8", "RApR7n_EM_CD8", "RAnR7n_EM_CD8")
    whldt$celltype <- factor(whldt$celltype, levels = x_celltype_levels)

    # Prepare new labels for facets panels
    cohort_labs <- c("HIV patients", "Healthy controls")
    names(cohort_labs) <- c('hiv', 'bcg')
    measurement_labs <- c("CCR5 positive cell propotion", "Mean of fluorescence intensity")
    names(measurement_labs) <- c('pc', 'gm')
    nlabs <- labeller(cohort = cohort_labs, measurement = measurement_labs)

    # Draw plot
    ggplot(data = whldt) +
      theme_bw() + 
      geom_point(aes_string(x = "celltype", y = "GeneSet", color = "-log10(adjP)"), size = 4) +
      facet_grid(cohort~measurement, labeller = nlabs) +
      xlab('Cell types') + ylab(cate) +
      scale_color_gradient(low = "grey", high = "red") +
      scale_y_discrete(limit = gene_sets, labels = y_ticks_labels) +
      theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust = 0.2, size = 10), axis.text.y = element_text(size = 10))

    # Save plot
    ggsave(paste0(cate, "-bubble_plot.png"), width = 11, height = 8)
    ggsave(paste0(cate, "-bubble_plot.pdf"), width = 11, height = 8)
    ggsave(paste0(cate, "-bubble_plot.svg"), width = 11, height = 8)
}

main()
