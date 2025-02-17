# Multi-Omics-Analysis-Workflow-with-MixOmics

# Data Processing and Analysis in R

This script processes and analyzes methylation, gene expression, and phenotype data using various R packages. It performs data cleaning, filtering, and normalization, followed by integration and visualization of the results.

## Required Libraries

The following R libraries are required for running this script:

```R
# Required libraries for data manipulation and analysis
library(data.table)     # For data table operations
library(tidyverse)      # For data wrangling and visualization
library(minfi)          # For methylation data analysis
library(zoo)            # For time series operations
library(limma)          # For linear modeling
library(lattice)        # For creating trellis graphics
library(ggplot2)        # For data visualization
library(reshape)        # For reshaping data
library(stringr)        # For string operations
library(RColorBrewer)   # For color palettes
library(mixOmics)       # For multivariate analysis (PLS-DA, etc.)
library(GeneNet)        # For network analysis
library(pROC)           # For ROC analysis
library(biomaRt)        # For gene annotation from Ensembl
library(org.Hs.eg.db)   # For human gene annotations
```

## Data Input

The following datasets are read into the script:

- `cleaned_train_meth.csv`: Training methylation data
- `cleaned_test_meth.csv`: Testing methylation data
- `cleaned_train_gene.csv`: Training gene expression data
- `cleaned_test_gene.csv`: Testing gene expression data
- `train_phenotype.csv`: Training phenotype data
- `test_phenotype.csv`: Testing phenotype data

The data are loaded using `read.csv()` and combined as follows:

```R
# Load methylation, gene expression, and phenotype data
train.meth = read.csv("cleaned_train_meth.csv")  # Read training methylation data
test.meth = read.csv("cleaned_test_meth.csv")    # Read testing methylation data

train.gene = read.csv("cleaned_train_gene.csv")  # Read training gene expression data
test.gene = read.csv("cleaned_test_gene.csv")    # Read testing gene expression data

train.phynotype = read.csv("train_phenotype.csv")  # Read training phenotype data
test.phynotype = read.csv("test_phenotype.csv")    # Read testing phenotype data

# Combine training and testing datasets
all.meth = cbind(train.meth, test.meth)
all.gene = cbind(train.gene, test.gene)
all.phynotype = rbind(train.phynotype, test.phynotype)
```

## Data Cleaning and Filtering

### Methylation Data

1. Removes the first column, sets row names, and omits missing values.
2. Filters out rows where all values are zero, removes rows with equal variance, and filters out rows with more than 90% zeros.

```R
# Remove first column and set row names for methylation data
all.meth <- all.meth[,-1]                       # Remove the first column
rownames(all.meth) <- all.meth[,1]               # Set row names to the first column
all.meth <- all.meth[,-1]                       # Remove the first column again
all.meth <- na.omit(all.meth)                    # Remove rows with missing values

# Filter out rows where all values are zero
rownames_all_meth <- rownames(all.meth)          # Store original row names
all.meth <- all.meth[rowSums(all.meth != 0, na.rm = TRUE) > 0, ]  # Remove rows with all zeros

# Filter rows with equal variance and more than 90% zeros
row_variance <- apply(all.meth, 1, var, na.rm = TRUE)  # Calculate variance for each row
all.meth_filtered <- all.meth[row_variance != row_variance[1], ]  # Remove rows with equal variance
numeric_data_meth <- all.meth_filtered[, sapply(all.meth_filtered, is.numeric)]  # Select numeric columns
zero_percentage_meth <- rowSums(numeric_data_meth == 0, na.rm = TRUE) / ncol(numeric_data_meth)  # Calculate zero percentage
all.meth_filtered <- all.meth_filtered[zero_percentage_meth <= 0.9, ]  # Filter rows with >90% zeros
```

### Gene Expression Data

1. Removes the first column, sets row names, and omits missing values.
2. Converts the data to numeric and filters out rows with all zeros, removes rows with equal variance, and filters out rows with more than 90% zeros.

```R
# Remove first column and set row names for gene expression data
all.gene <- all.gene[,-1]                       # Remove the first column
rownames(all.gene) <- all.gene$V1               # Set row names to the first column
all.gene <- all.gene[,-1]                       # Remove the first column again
all.gene <- na.omit(all.gene)                    # Remove rows with missing values

# Convert data to numeric and round values
rownames_all_gene <- rownames(all.gene)          # Store original row names
all.gene <- as.data.frame(sapply(all.gene, as.numeric))  # Convert to numeric
all.gene <- round(all.gene)                      # Round the values
rownames(all.gene) <- rownames_all_gene          # Set the row names again

# Filter out rows where all values are zero
all.gene <- all.gene[rowSums(all.gene != 0, na.rm = TRUE) > 0, ]  # Remove rows with all zeros

# Filter rows with equal variance and more than 90% zeros
row_variance <- apply(all.gene, 1, var, na.rm = TRUE)  # Calculate variance for each row
all.gene_filtered <- all.gene[row_variance != row_variance[1], ]  # Remove rows with equal variance
numeric_data <- all.gene[, sapply(all.gene, is.numeric)]  # Select numeric columns
zero_percentage <- rowSums(all.gene_filtered == 0, na.rm = TRUE) / ncol(all.gene_filtered)  # Calculate zero percentage
all.gene_filtered <- all.gene_filtered[zero_percentage <= 0.9, ]  # Filter rows with >90% zeros
```

### Phenotype Data

Removes the `sample` column and sets row names:

```R
# Remove 'sample' column and set row names for phenotype data
rownames(all.phynotype) <- all.phynotype$sample  # Set row names to 'sample' column
all.phynotype <- all.phynotype[, -which(names(all.phynotype) == "sample")]  # Remove 'sample' column
```

### Aligning Datasets by Sample ID

Aligns the rows of both the gene and methylation data by matching sample IDs:

```R
# Align methylation and gene expression data by matching sample IDs
sample_ids_gene <- colnames(all.gene_filtered)  # Get gene expression sample IDs
sample_ids_meth <- colnames(all.meth_filtered)  # Get methylation sample IDs
common_samples <- intersect(colnames(all.gene_filtered), colnames(all.meth_filtered))  # Find common samples
length(common_samples)  # Print out the number of common samples

# Subset the data to include only common samples
all.gene_filtered <- all.gene_filtered[, common_samples, drop = FALSE]
all.meth_filtered <- all.meth_filtered[, common_samples, drop = FALSE]

# Check if the number of rows is greater than zero
nrow(all.gene_filtered)  # Should be > 0
nrow(all.meth_filtered)  # Should be > 0
```

### Gene ID Conversion

Maps Ensembl gene IDs to gene symbols and resolves duplicates:

```R
# Convert Ensembl IDs to gene symbols
ensembl_ids_no_version <- sub("\\..*", "", rownames(all.gene_filtered))  # Remove version from Ensembl IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids_no_version, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")  # Map to gene symbols

# Handle duplicated gene symbols
duplicated_genes <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)  # Identify duplicates
gene_symbols[duplicated_genes] <- paste0(gene_symbols[duplicated_genes], "_", which(duplicated_genes))  # Resolve duplicates
length(unique(gene_symbols)) == length(gene_symbols)  # Check if there are no duplicates
rownames(all.gene_filtered) <- gene_symbols  # Update row names with gene symbols
```
---

# Data Analysis and Visualization

This section of the analysis focuses on various statistical and machine learning methods for evaluating, integrating, and visualizing gene expression and methylation data. It includes performance evaluation, Principal Component Analysis (PCA), Partial Least Squares (PLS) analysis, Sparse PLS, PLS Discriminant Analysis (PLS-DA), and DIABLO integration.

## Performance Evaluation with Cross-Validation

To evaluate the model's performance, we perform cross-validation using the `perf()` function with 5-fold validation and 10 repeats. This process helps assess how well the PLS-DA model generalizes to unseen data.

```R
# Set seed for reproducibility
set.seed(30)

# Performance evaluation using cross-validation
perf_plsda <- perf(plsda_model, validation = 'Mfold', folds = 5, 
                   progressBar = FALSE, nrepeat = 10)

# Plot performance results with standard deviation and legend
plot(perf_plsda, sd = TRUE, legend.position = 'horizontal')
```

---

## Principal Component Analysis (PCA)

PCA is performed to reduce the dimensionality of both the gene expression and methylation data, while capturing the most significant variation. The first few principal components are visualized for both datasets.

### PCA on Gene Expression

```R
# Perform PCA on gene expression data with 4 components
pca.gene <- pca(t_all.gene_filtered, ncomp = 4)

# Plot the PCA result for gene expression
plotIndiv(pca.gene, group = groups, legend = TRUE, title = "Gene Expression PCA")
```

### PCA on Methylation

```R
# Perform PCA on methylation data with 4 components
pca.meth <- pca(t_all.meth_filtered, ncomp = 4)

# Plot the PCA result for methylation data
plotIndiv(pca.meth, group = groups, legend = TRUE, title = "Methylation Expression PCA")
```

---

## Partial Least Squares (PLS) Integration

PLS is used to model the relationship between gene expression and methylation data, and we visualize the integration of these two data types.

```R
# PLS integration between gene and methylation data with 4 components
pls.result_1 <- pls(X = t_all.gene_filtered, Y = t_all.meth_filtered, ncomp = 4)

# Plot the integrated PLS result
plotIndiv(pls.result_1, group = groups, legend = TRUE, title = "PLS Integration")
```

---

## Sparse Partial Least Squares (sPLS)

Sparse PLS is employed to identify the most relevant variables (features) in the gene expression data. In this case, we keep 50 variables for both X (gene expression) and Y (gene expression).

```R
# Perform Sparse PLS with 4 components and 50 variables for both X and Y
spls.result <- spls(X = t_all.gene_filtered, Y = t_all.gene_filtered, ncomp = 4, keepX = c(50, 50), keepY = c(50, 50))

# Plot the sPLS correlation circle
plotVar(spls.result, cutoff = 0.8, title = "sPLS Correlation Circle", cex = c(6, 6))
```

---

## PLS-DA on Gene Expression

PLS-DA is performed on the gene expression data to classify samples based on phenotype groups. This method helps in identifying key features that discriminate between different sample types.

```R
# Perform PLS-DA on gene expression data with 4 components
plsda.result <- plsda(X = t_all.gene_filtered, Y = groups, ncomp = 4)

# Plot the PLS-DA result
plotIndiv(plsda.result, group = groups, legend = TRUE, title = "PLS-DA on Gene Expression")
```

---

## DIABLO Integration (for Multi-Omics Data)

DIABLO (Data Integration Analysis for Biomarker Discovery using Latent components) integrates multiple omics data, such as gene expression and methylation, into a unified analysis. The results are visualized to identify significant components from both datasets.

```R
# DIABLO integration of gene and methylation data with 4 components
diablo.result <- block.plsda(X = list(gene = t_all.gene_filtered, meth = t_all.meth_filtered), Y = groups, ncomp = 4)

# Plot the DIABLO integration result
plotIndiv(diablo.result, legend = TRUE, title = "DIABLO Integration")

# Plot the DIABLO variable contributions
plotVar(diablo.result, cutoff = 0.85, title = "DIABLO Variable Plot", cex = c(6, 6))
```

---

## Thanks!
---
