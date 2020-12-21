#!/usr/bin/env Rscript
library('ggplot2')
library('optparse')
library('data.table')
library('reshape')

getargs <- function() {
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), dest="input", help="Input file.")
  parser <- add_option(parser, c("-o", "--output"), dest="output", help="Output file.")

  return(parse_args(parser))
}

my_cast <- function(..., rnm_col=NULL) {
  res <- reshape::cast(...)

  if (!is.null(rnm_col)) {
    rownames(res) <- res[, rnm_col]
    res <- res[, colnames(res) != rnm_col]
  }
  res
}

main <- function() {
  opts <- getargs()
  input <- opts$input
  output <- opts$output

  if (!file.exists(input)) { stop(paste0("Not found ", input)) }
  tbl <- fread(input, data.table=FALSE)

  cttbl <- sort(table(tbl[, "CellType1"]))
  ct <- names(cttbl)

  tbl <- my_cast(tbl, CellType1 ~ CellType2, value = "PP4", rnm_col="CellType1")

  if (nrow(tbl) != length(ct)) {
    tbl[ct[!ct %in% rownames(tbl)], ] <- NA
  }
  print(dim(tbl))

  if (ncol(tbl) != length(ct)) {
    tbl[, ct[!ct %in% colnames(tbl)]] <- NA
  }
  print(dim(tbl))

  print(tbl[ct, ct])

  # ct_sub <- c(ct[ct %in% names(tbl)])
  # print(tbl[ct, ct_sub])
}

main()
