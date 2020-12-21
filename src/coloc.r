#!/usr/bin/env Rscript

library(coloc)
library(stringi)
library(optparse)
library(data.table)
library(zzhRmisc)

# CHR BP SNP A1 A2 N P P(R) OR OR(R) Q I

caw <- function(dataset1, dataset2, maf=NULL) {
#A wrapper of coloc.abf()
  if (is.null(maf))
    return(coloc.abf(dataset1=dataset1, dataset2=dataset2))

  return(coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=maf))
}


mkds <- function(dataset, type='quant', pval_col=NULL, beta_col=NULL,
                 var_col=NULL, y_col=NULL, nsample=NULL, maf=NULL) {
#Make dataset for coloc_abf_wapper()
    if (!is.null(pval_col))
        return(list(pvalues=dataset[, pval_col], N=nsample, type=type, MAF=maf))

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
  ssdf <- fread(ssfp, data.table=FALSE, header=TRUE)
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


# Get options
getargs <- function() {
  parser <- OptionParser(description='A wrapper of coloc.abf for p-values.')
  
  parser <- add_option(parser, c('-a', '--ctrl'), type='character', dest='ctrl_ss', metavar='CTRL-SUM-STAT', help='The summary statistics of control cohort.')
  parser <- add_option(parser, c('-b', '--case'), type='character', dest='case_ss', metavar='CASE-SUM-STAT', help='The summary statistics of case cohort.')
  parser <- add_option(parser, c('-d', '--ctrl-maf'), type='character', dest='ctrl_maf_fl', metavar='MAF-FILE', help='Dataset of minor allele afrequency for summary statistics -a.')
  parser <- add_option(parser, c('-e', '--case-maf'), type='character', dest='case_maf_fl', metavar='MAF-FILE', help='Dataset of minor allele afrequency for summary statistics -b.')

  parser <- add_option(parser, c('-A', '--ctrl-cols'), dest='ctrl_cols', metavar='SNP,PVL,CHR,POS', default='snps,pvalue,SequenceName,Position', help='The column ids <SNP,PVL,CHR,POS> for ctrl summary statistics. Default: %default')
  parser <- add_option(parser, c('-B', '--case-cols'), dest='case_cols', metavar='SNP,PVL,CHR,POS', default='snps,pvalue,SequenceName,Position', help='The column ids <SNP,PVL,CHR,POS> for case summary statistics. Default: %default')
  parser <- add_option(parser, c('-D', '--ctrl-maf-cols'), type='character', dest='ctrl_maf_cols', metavar='SNP,MAF,CHR,POS', default='snps,MAF,SequenceName,Position', help='The column ids <SNP,MAF,CHR,POS> for MAF files. Default: %default')
  parser <- add_option(parser, c('-E', '--case-maf-cols'), type='character', dest='case_maf_cols', metavar='SNP,MAF,CHR,POS', default='snps,MAF,SequenceName,Position', help='The column ids <SNP,MAF,CHR,POS> for MAF files. Default: %default')

  parser <- add_option(parser, c('-c', '--chrom'), type='character', dest='chrom', metavar='CHR', default='3', help='The target chromosome.')
  parser <- add_option(parser, c('--max-pos'), type='integer', dest='max_pos', metavar='MAX-POS', default=NULL, help='The end position of the region.')
  parser <- add_option(parser, c('--min-pos'), type='integer', dest='min_pos', metavar='MIN-POS', default=NULL, help='The start position of the region.')

  parser <- add_option(parser, c('-z', '--smpl-size'), type='character', dest='smpl_size', metavar='SAMPLE-SIZE,SMAPLE-SIZE', help='The sample size.')
  
  return(parse_args(parser))
}


main <- function() {
  opts <- getargs()
  
  if (is.null(opts$ctrl_ss)) stop('-a/--ctrl is required!')
  if (is.null(opts$case_ss)) stop('-b/--case is required!')
  if (is.null(opts$case_maf_fl)) stop('-d/--ctrl-maf is required!')
  if (is.null(opts$case_maf_fl)) stop('-e/--case-maf is required!')
  if (is.null(opts$smpl_size)) stop('-z/--smpl-size is required!')

  chrom <- opts$chrom
  ctrl_maf_fl <- opts$ctrl_maf_fl
  case_maf_fl <- opts$case_maf_fl

  c(ct_snp, ct_pvl, ct_chr, ct_pos) %=% stri_split(opts$ctrl_cols, regex=',')[[1]]
  c(cs_snp, cs_pvl, cs_chr, cs_pos) %=% stri_split(opts$case_cols, regex=',')[[1]]
  c(ctmf_snp, ctmf_maf, ctmf_chr, ctmf_pos) %=% stri_split(opts$ctrl_maf_cols, regex=',')[[1]]
  c(csmf_snp, csmf_maf, csmf_chr, csmf_pos) %=% stri_split(opts$case_maf_cols, regex=',')[[1]]
  c(ct_sm_size, cs_sm_size) %=% as.integer(stri_split(opts$smpl_size, regex=',')[[1]])

  ctrl_ss <- lddt(opts$ctrl_ss, ct_snp, ct_chr, chrom)
  case_ss <- lddt(opts$case_ss, cs_snp, cs_chr, chrom)
  ctrl_maf <- lddt(opts$ctrl_maf_fl, ctmf_snp, ctmf_chr, chrom)
  case_maf <- lddt(opts$case_maf_fl, csmf_snp, csmf_chr, chrom)

  ss_cand_snps <- intersect(rownames(ctrl_ss), rownames(case_ss))
  mf_cand_snps <- intersect(rownames(ctrl_maf), rownames(case_maf))
  cand_snps <- intersect(ss_cand_snps, mf_cand_snps)

  ctrldt <- mkds(ctrl_ss[cand_snps, ], pval_col=ct_pvl, nsample=ct_sm_size,
                 maf=ctrl_maf[cand_snps, ctmf_maf])
  casedt <- mkds(case_ss[cand_snps, ], pval_col=cs_pvl, nsample=cs_sm_size,
                 maf=case_maf[cand_snps, csmf_maf])

  res <- caw(ctrldt, casedt)
}

main()
