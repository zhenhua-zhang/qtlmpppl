#!/usr/bin/env Rscript
suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(TwoSampleMR))


get_opt <- function() {
    parser <- OptionParser(description = "A wrapper of TwoSampleMR.")
    parser <- add_option(parser, c("-e", "--exposure"), dest = "exp_file", help = "The exposure file. Summary statistics from MatrixEQTL.")
    parser <- add_option(parser, c("-o", "--outcome"), dest = "out_file", help = "The outcome file. Summary satatistics from MatrixEQTL.")
    parser <- add_option(parser, c("-a", "--annot-info"), dest = "annot_info", help = "Annotation information in VCF format (could be gzip compressed).")
    
    return(parse_args(parser))
}


se_by_pb <- function(pv, bt, ss, np) {
    # Calculate SE from p-value and beta (conefficient)
    return(abs(bt/qt(pv/2, df = ss - np)))
}


load_annot <- function(fpath) {
    # Load annotation data in VCF format
    vcf <- read.vcfR(fpath, verbose = FALSE)

    return(vcf)
}


load_ssdat <- function() {
    # Load summary statistics data in csv or tsv.
}


annot_ss_dat <- function(ss_dat, annot) {
}


main <- function() {
    opts <- get_opt()
    exp_file <- opts$exp_file
    out_file <- opts$out_file
    annot_info <- opts$annot_info

    load_annot(annot_info)

    # exp_dat <- read_exposure_data(exp_file)
    # out_dat <- read_outcome_data(out_file)

    # har_dat <- harmonise_data(exp_dat, out_dat)
    # mr_res <- mr(har_dat)
}

main()
