# QTL mapping pipeline for CCR5 expression in 200HIV project

## Package tree

```

```


# Genotype

Questions:
    1. Is it necessary to do the imputation by combining BCG and HIV cohort?

## Genotyping

## Imputation

## Quality control

# Phenotype

## Preprocessing
1. Convert data from excel into `csv` or `tsv`
2. Check distribution and outliers
    - Check distribution by histogram, transform the data by `log2` or `log10`, `sqrt` if necessary
    - Mask outlier measurements as NA per trait

## Correlation analysis
1. Using Partial Correlation Analysis to generate a heatmap while adjusting for other co-variables. Thanks Louise's inspiration.


# QTL mapping analysis
1. `MatrixEQTL`
2. Manhattan plot. E.g. `circos`
3. QQ-plot, and inflation λ

# Function Annotation

# Replication of

## CCR5AS
**Description** : CCR5 antisense RNA
**BioType**     : LncRNA
**Associations**
    rs6441975-?
    rs62246129-T
    rs3087253-?
    rs2373226-T
    rs2213290-T
    rs1001007-?
    rs1800024-T  # Not in the genotypes of 200HIV cohort
    rs62625034-G # Not in the genotypes of 200HIV cohort

## CCR5
**Description** : C-C motif chemokine receptor 5(gene/pseudogene)
**Biotype**     : protein_coding
**Associations**
    rs333(indel)

## HIV-1 virus load
**Reference**: Polymorphisms of large effect explain the majority of the host genetic contribution of variation of HIV-1 virus load. PNAS, 2015
**Associations**
    rs1015164 (full samples)
    rs4317138 (CCR5Δ32 and Hap-P1 carriage)

## Epigenetic

## Reported SNPs



# Note
1. In 300BCG cohort, the measurement of GM of Naive_CD8 is not normal respect
   to normal distribution even converted by log10.



# Results

Date: 2020 Sep 07 13:20:22
[x] 1. QTLs for each association analysis
[ ] 2. Check if the associations are the same locus
    [x] - Coloc
    [ ] - PLINK clump
3.
