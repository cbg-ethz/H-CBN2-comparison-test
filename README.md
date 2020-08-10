# H-CBN2-comparison-test

Distance-based test used to compare mutational networks estimated using the H-CBN2 method

## Dependencies

- Snakemake

The following R packages are required:
- mccbn
- optparser

## Usage

### Input
Two genotype data sets saved as RDS files. Each RDS file contains a matrix with N rows (N: number of genotypes) and p columns (p: number of mutations).

### Output
The output includes the Jaccard distance between estimated posets, the corresponding p-value, and the empirical distribution of Jaccard distances obtained through a permutation test.

### Example
```
snakemake -p --cores 1
```
