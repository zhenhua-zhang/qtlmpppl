#!/usr/bin/env Rscript
# A script to analyze KEGG profile for genes mapped by CCR5 QTLs
# R: 3.5.2
# Packages:
#   Bioconductor: 3.8
#   clusterProfiler: 
#   enrichPlot:
# Don't forget to cite:
#   Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R
#   package for comparing biological themes among gene clusters. OMICS: A
#   Journal of Integrative Biology. 2012, 16(5):284 -287.

print("The analysis looks less useful as there is a few genes (65) to analyze.")

setwd("~/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/mapped_genes")

library(clusterProfiler)
library(enrichplot)

# Ensure the organism is Magnaporthe oryzae
# search_kegg_organism("mgr", by = "kegg_code")

ko_list <- read.table("mapped_genes.txt", header = 1)

genes = ko_list$entrez
organ = "hsa"
show_n = 20
kpe <- enrichKEGG(gene = genes, organism = organ, pvalueCutoff = 0.05, pAdjustMethod = "fdr")

write.csv(kpe, "./clusterProfiler/kegg_pathway.enrichment.csv")

# barplot(kpe, showCategory = show_n)
# heatplot(kpe, showCategory = show_n)
pdf("./clusterProfiler/kegg_pathway.upsetplot.pdf", width = 15, height = 10)
upsetplot(kpe, n = show_n)
dev.off()

pdf("./clusterProfiler/kegg_pathway.dotplot.pdf")
dotplot(kpe, showCategory = show_n)
dev.off()

pdf("./clusterProfiler/kegg_pathway.enrich_map.pdf")
emapplot(kpe, showCategory = show_n)
dev.off()
