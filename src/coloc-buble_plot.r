#!/usr/bin/env Rscript

library(dplyr)
library(stringi)
library(ggplot2)
library(gridExtra)
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
        bsfn <- stri_replace(fn, "", fixed = "-per_cell_type.txt")
        tmpdt <- fread(fn, verbose = FALSE)
        tmpdt[, `:=`(comparison = bsfn)]
        
        if (is.null(whldt)) { whldt <- tmpdt }
        else { whldt <- rbind(whldt, tmpdt) }
    }
    
    return(whldt)
}


main <- function() {
    # Load data and select category
    chosen_cp <- c("CD45", "L", "CD4", "mTreg", "TEM", "RApR7n", "RAnR7n", "CD8", "CM_CD8", "EM_CD8", "RApR7n_EM_CD8", "RAnR7n_EM_CD8")
    whldt <- lddt("*-per_cell_type.txt") %>%
      filter(Chrom == 3, CellType %in% chosen_cp)

    whldt$CellType <- factor(whldt$CellType, levels = chosen_cp)

    # Prepare new labels for facets panels
    nlabs <- labeller(Cohort = c("hiv" = "HIV patients", "bcg" = "Healthy controls"))

    # Draw plot
    p_ccr5Exp_viralLoad <- whldt %>%
      filter(comparison == "ccr5-vs-viral_load", CellType %in% chosen_cp) %>%
      rowwise() %>%
      mutate(Cohort = stri_split_fixed(Measurement, pattern = "-", n = 2)[[1]][1],
             Measurement = stri_split_fixed(Measurement, pattern = "-", n = 2)[[1]][2],
             "PP4>=0.8" = PP4 >= 0.8) %>%
      as.data.table() %>%
      ggplot() +
        theme_bw() + 
        geom_point(aes_string(x = "Measurement", y = "CellType", color = "PP4", shape = "PP4>=0.8"), size = 5) +
        facet_grid(cols=vars(Cohort), labeller = nlabs) +
        ylab(element_blank()) +
        xlab(element_blank()) +
        scale_color_gradient(low = "grey", high = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust = 0.2, size = 10),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    # Draw plot
    p_cohort <- whldt %>%
      filter(comparison == "comparison_cohort") %>%
      rowwise() %>%
      mutate("PP4>=0.8" = PP4 >= 0.8, "cmp_type" = "Cohort") %>%
      as.data.table() %>%
      ggplot() +
        theme_bw() +
        geom_point(aes_string(x = "Measurement", y = "CellType", color = "PP4", shape = "PP4>=0.8"), size = 5) +
        facet_grid(cols=vars(cmp_type)) +
        ylab('Cell types') +
        xlab(element_blank()) +
        scale_color_gradient(low = "grey", high = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = .95, vjust = 0.2, size = 10),
              legend.position = "none")


    # Draw plot
    p_measurement <- whldt %>%
      filter(comparison == "comparison_measurements") %>%
      rowwise() %>%
      mutate("PP4>=0.8" = PP4 >= 0.8, "cmp_type" = "Measurement") %>%
      as.data.table() %>%
      ggplot() +
        theme_bw() +
        geom_point(aes_string(x = "Measurement", y = "CellType", color = "PP4", shape = "PP4>=0.8"), size = 5) +
        facet_grid(cols=vars(cmp_type)) +
        ylab(element_blank()) +
        xlab(element_blank()) +
        scale_color_gradient(low = "grey", high = "red") +
        scale_x_discrete(labels = c("hiv" = "HIV patients", "bcg" = "Healthy Controls")) +
        theme(axis.text.x = element_text(angle = 90, hjust = .9, vjust = 0.2, size = 10),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "none")

    p <- grid.arrange(p_cohort, p_measurement, p_ccr5Exp_viralLoad, nrow = 1, widths = c(1.7, 1, 2.5))

    # Save plot
    ggsave("coloc-bubble_plot.pdf", p, width = 8, height = 6)
    ggsave("coloc-bubble_plot.svg", p, width = 8, height = 6)
    ggsave("coloc-bubble_plot.png", p, width = 8, height = 6)
}

main()
