#!/usr/bin/env Rscript

library(coloc)
library(stringi)
library(optparse)
library(data.table)
library(zzhRmisc)

# CHR BP SNP A1 A2 N P P(R) OR OR(R) Q I


# Get options
getargs <- function() {
  parser <- OptionParser(description = 'A wrapper of coloc.abf for p-values.')
  
  parser <- add_option(parser, c('-a', '--ctrl'), type = 'character', dest = 'ctrl_ss', metavar = 'CTRL-SUM-STAT', help = 'The summary statistics of control cohort.')
  parser <- add_option(parser, c('-b', '--case'), type = 'character', dest = 'case_ss', metavar = 'CASE-SUM-STAT', help = 'The summary statistics of case cohort.')
  parser <- add_option(parser, c('-m', '--mafp'), type = 'character', dest = 'maf_fl', metavar = 'MAF-FILE', help = 'Dataset of minor allele afrequency.')

  parser <- add_option(parser, c('-A', '--ctrl-cols'), dest = 'ctrl_cols', metavar = 'SNP,PVL,CHR,POS', default = 'snps,pvalue,SequenceName,Position', help = 'The column ids <SNP,PVL,CHR,POS> for ctrl summary statistics. Default: %default')
  parser <- add_option(parser, c('-B', '--case-cols'), dest = 'case_cols', metavar = 'SNP,PVL,CHR,POS', default = 'snps,pvalue,SequenceName,Position', help = 'The column ids <SNP,PVL,CHR,POS> for case summary statistics. Default: %default')
  parser <- add_option(parser, c('-M', '--maff-cols'), dest = 'maff_cols', metavar = 'SNP,MAF,CHR,POS', default = 'snps,MAF,SequenceName,Position', help = 'The column ids <SNP,MAF,CHR,POS> for MAF files. Default: %default')

  parser <- add_option(parser, c('-c', '--chrom'), type = 'character', dest = 'chrom', metavar = 'CHR', default = '3', help = 'The target chromosome.')
  parser <- add_option(parser, c('--pos-max'), type = 'integer', dest = 'max_pos', metavar = 'MAX-POS', default = NULL, help = 'The end position of the region.')
  parser <- add_option(parser, c('--pos-min'), type = 'integer', dest = 'min_pos', metavar = 'MIN-POS', default = NULL, help = 'The start position of the region.')

  parser <- add_option(parser, c('-z', '--sm-size'), type = 'character', dest = 'sm_size', metavar = 'SAMPLE-SIZE,SMAPLE-SIZE', help = 'The sample size.')
  
  return(parser)
}


caw <- function(dataset1, dataset2, maf = NULL) {
#A wrapper of coloc.abf()
  if (is.null(maf))
    return(coloc.abf(dataset1 = dataset1, dataset2 = dataset2))

  return(coloc.abf(dataset1 = dataset1, dataset2 = dataset2, MAF=maf))
}


mkds <- function(dataset, type='quant', pval_col=NULL, beta_col=NULL,
                 var_col=NULL, y_col=NULL, nsample=NULL) {
#Make dataset for coloc_abf_wapper()
    if (!is.null(pval_col))
        return(list(pvalues=dataset[, pval_col], N=nsample, type=type))

    if ((!is.null(beta_col)) && (!is.null(var_col)) && (!is.null(y_col)))
        return(list(beta=dataset[, beta_col], varbeta=dataset[, var_col],
                    N=nsample, sdY=sd(dataset[, y_col]), type=type))
    else if (is.null(y_col))  # FIXME: not a correct way to do ABF analysis
        return(NULL)
        #return(list(beta=dataset[, beta_col], varbeta=dataset[, var_col],
                    #N=nsample, type=type))

    return(NULL)
}


lddt <- function(ssfp, snp_col, chr_col=NULL, chr=NULL, pos_col=NULL,
                 pos=c(NULL, NULL)) {
#Load input file into a data.table
  ssdf <- fread(ssfp, data.table = FALSE, header = TRUE)
  ssdf <- ssdf[which(!ssdf[, snp_col] %in% c('.')), ]
  rownames(ssdf) <- ssdf[, snp_col]

  if (is.null(pos_col)) {
    pos_min_cond <- TRUE
    pos_max_cond <- TRUE
  } else {
    min_pos <- pos[1]
    max_pos <- pos[2]

    if (is.null(min_pos))
        pos_min_cond <- TRUE
    else
        pos_min_cond <- ssdf[, pos_col] >= min_pos

    if (is.null(max_pos))
        pos_max_cond <- TRUE
    else
        pos_max_cond <- ssdf[, pos_col] <= max_pos
  }

  ssdf[, chr_col] <- as.character(ssdf[, chr_col])
  if (is.null(chr_col)) {
    chr_cond <- TRUE
  } else {
    if (is.null(chr_col))
        chr_cond <- TRUE
    else
        chr_cond <- ssdf[, chr_col] == chr
  }

  return(ssdf[which(chr_cond & pos_max_cond & pos_min_cond), ])
}


main <- function() {
  opts <- parse_args(getargs())
  
  ctrl_ss <- opts$ctrl_ss
  if (is.null(ctrl_ss))
    stop('-a/--ctrl is required!')

  case_ss <- opts$case_ss
  if (is.null(case_ss))
    stop('-b/--case is required!')

  chrom <- opts$chrom
  maf_fl <- opts$maf_fl

  c(ct_snp_col, ct_pvl_col, ct_chr_col, ct_pos_col) %=% stri_split(opts$ctrl_cols, regex=',')[[1]]
  c(cs_snp_col, cs_pvl_col, cs_chr_col, cs_pos_col) %=% stri_split(opts$case_cols, regex=',')[[1]]
  c(mf_snp_col, mf_maf_col, mf_chr_col, mf_pos_col) %=% stri_split(opts$maff_cols, regex=',')[[1]]
  c(ct_sm_size, cs_sm_size) %=% as.integer(stri_split(opts$sm_size, regex=',')[[1]])

  ctrldf <- lddt(ctrl_ss, ct_snp_col, ct_chr_col, chrom)
  casedf <- lddt(case_ss, cs_snp_col, cs_chr_col, chrom)

  cand_snps <- intersect(rownames(ctrldf), rownames(casedf))

  maf <- NULL
  if (!is.null(maf_fl)) {
    maf_df <- lddt(maf_fl, mf_snp_col, mf_chr_col, chrom)
    cand_snps <- intersect(cand_snps, rownames(maf_df))
    maf <- maf_df[cand_snps, mf_maf_col]
    ctrldt <- mkds(ctrldf[cand_snps, ], pval_col=ct_pvl_col, nsample=ct_sm_size)
    casedt <- mkds(casedf[cand_snps, ], pval_col=cs_pvl_col, nsample=cs_sm_size)
  } 

  res <- caw(ctrldt, casedt, maf)
}

main()
