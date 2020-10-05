#!/usr/bin/env Rscript
library(LDheatmap)


# Working dir
setwd("~/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/significant_snps/")


snps_pos <- unique(read.csv('./significant_snps.txt', header = 1)[, c("snps", "Position")])
rownames(snps_pos) <- snps_pos[, "snps"]

lddtfm <- read.table("./chr3-nodup-ld_matrix/significant_snps-chr3-nodup-ldmatrix.txt", header = 1, row.names = 1)
ld_snps <- rownames(lddtfm)

snps_pos_ld <- snps_pos[ld_snps, "Position"]

rgb.palette <- colorRampPalette(rev(c("white", "grey", "red")), space = "rgb")
llh <- LDheatmap(as.matrix(lddtfm), genetic.distances = snps_pos_ld, color = rgb.palette(18), flip = T)
llh <- LDheatmap.addGenes(llh, chromosome = "chr3", genome = "hg18", non_coding = TRUE, splice_variants = FALSE)

pdf("./chr3-nodup-ld_matrix/ld_matrix-heatmap-gene.pdf")
llh <- LDheatmap.addScatterplot(llh, snps_pos_ld)
dev.off()