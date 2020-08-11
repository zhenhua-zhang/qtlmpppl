# QTL mapping pipeline for CCR5 expression in 200HIV project

## Package tree

```

~/Documents/projects/wp_hiv_ccr5/scripts
├── collect_and_report.py  [Collecting outputs of qtl_mapping.r and generate a report]
├── preprocessing.py       [Preprocessing input file for qtl_mapping.r]
├── qtl_mapping.r          [Do QTL mapping exploiting MatrixEQTL package]
├── README.md              [This readme file]
└── run.sh                 [Execute all steps in one script.]

0 directories, 5 files
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

|------------|------------------------|----------------|----------------------|-------------------|-------------------------------------------|---------|--------------------|------------------------------------|------------------------|------------------------------------|---------|-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| SNP ID     | Position               | Trait (200HIV) | p-value (200HIV)     | Beta (200HIV)     | Disease or trait (public)                 | p-value | Odds ratio or beta | Sample size (public)               | Reported Gene (public) | Strongest SNP-risk allele (public) | Alleles | Risk Allele Frequency | Publication                                                                                                                                                              |
|------------|------------------------|----------------|----------------------|-------------------|-------------------------------------------|---------|--------------------|------------------------------------|------------------------|------------------------------------|---------|-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| rs79815064 | chr3:46277577-46277577 | RAnR7n_EM_CD8  | 1.84958941322281e-10 | 0.627370236574362 | Macrophage inflammatory protein 1b levels | 3E-126  | 0.5505             | 8,243 Finnish ancestry individuals | Not reported           | rs79815064-G                       | A/G     | Not reported          | Ahola-Olli AV et al. Genome-wide Association Study Identifies 27 Loci Influencing Concentrations of Circulating Cytokines and Growth Factors. Am J Hum Genet. 2016-12-13 |
| rs4317138  |                        |                |                      |                   | HIV-1 virus load                          |         |                    |                                    |                        |                                    |         |                       | Polymorphisms of large effect explain the majority of the host genetic contribution of variation of HIV-1 virus load                                                     |
| rs1015164  |                        |                |                      |                   | HIV-1 virus load                          |         |                    |                                    |                        |                                    |         |                       | Polymorphisms of large effect explain the majority of the host genetic contribution of variation of HIV-1 virus load                                                     |
|------------|------------------------|----------------|----------------------|-------------------|-------------------------------------------|---------|--------------------|------------------------------------|------------------------|------------------------------------|---------|-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|



# Note
1. In 300BCG cohort, the measurement of GM of Naive_CD8 is not normal respect
   to normal distribution even converted by log10.
