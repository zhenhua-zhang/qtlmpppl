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
