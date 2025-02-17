# Load required libraries
library(data.table)  # For fast data manipulation
library(tidyverse)   # For data manipulation and visualization
library(minfi)       # For methylation data analysis
library(zoo)         # For time series and other irregular data handling
library(limma)       # For differential expression analysis
library(lattice)     # For data visualization using trellis
library(ggplot2)     # For data visualization
library(reshape)     # For reshaping data
library(stringr)     # For string manipulation
library(RColorBrewer) # For color palettes
library(mixOmics)    # For multivariate data analysis
library(GeneNet)     # For network analysis of gene expression data
library(pROC)        # For ROC analysis
library(biomaRt)     # For accessing genomic data
library(org.Hs.eg.db) # For gene annotation

# Read input files
train.meth=read.csv("cleaned_train_meth.csv")   # Train methylation data
test.meth=read.csv("cleaned_test_meth.csv")     # Test methylation data

train.gene=read.csv("cleaned_train_gene.csv")   # Train gene expression data
test.gene=read.csv("cleaned_test_gene.csv")     # Test gene expression data

train.phynotype=read.csv("train_phenotype.csv") # Train phenotype data
test.phynotype=read.csv("test_phenotype.csv")   # Test phenotype data

# Combine train and test datasets
all.meth=cbind(train.meth,test.meth)            # Combine methylation data
all.gene=cbind(train.gene,test.gene)            # Combine gene expression data
all.phynotype=rbind(train.phynotype,test.phynotype) # Combine phenotype data

###########################################################
# Data cleaning for methylation data
all.meth <- all.meth[,-1]  # Remove first column (row names)
rownames(all.meth) <- all.meth[,1]  # Set row names from the first column
all.meth <- all.meth[,-1]  # Remove the first column again
all.meth <- na.omit(all.meth)  # Remove rows with missing values

# Save row names for later use
rownames_all_meth <- rownames(all.meth)

# Remove rows with all zero values
all.meth <- all.meth[rowSums(all.meth != 0, na.rm = TRUE) > 0, ]

# Calculate row variance and filter out rows with equal variance
row_variance <- apply(all.meth, 1, var,na.rm = TRUE)
all.meth_filtered <- all.meth[row_variance != row_variance[1], ]

# Filter rows with more than 90% zeros
numeric_data_meth <- all.meth_filtered[, sapply(all.meth_filtered, is.numeric)]
zero_percentage_meth <- rowSums(numeric_data_meth == 0, na.rm = TRUE) / ncol(numeric_data_meth)
all.meth_filtered <- all.meth_filtered[zero_percentage_meth <= 0.9, ]

###################################################################################

# Data cleaning for gene expression data
all.gene<-all.gene[,-1]  # Remove first column (row names)
rownames(all.gene)<-all.gene$V1  # Set row names from first column
all.gene<-all.gene[,-1]  # Remove the first column again
all.gene <- na.omit(all.gene)  # Remove rows with missing values

# Save row names for later use
rownames_all_gene <- rownames(all.gene)
all.gene <- as.data.frame(sapply(all.gene, as.numeric))  # Convert all columns to numeric

# Round the data
all.gene <- round(all.gene)
rownames(all.gene) <- rownames_all_gene

# Remove rows with all zero values
all.gene <- all.gene[rowSums(all.gene != 0, na.rm = TRUE) > 0, ]

# Calculate row variance and filter out rows with equal variance
row_variance <- apply(all.gene, 1, var, na.rm = TRUE)
all.gene_filtered <- all.gene[row_variance != row_variance[1], ]

# Filter rows with more than 90% zeros
numeric_data <- all.gene[, sapply(all.gene, is.numeric)]
zero_percentage <- rowSums(all.gene_filtered == 0, na.rm = TRUE) / ncol(all.gene_filtered)
all.gene_filtered <- all.gene_filtered[zero_percentage <= 0.9, ]

################################################################################
# Clean phenotype data
rownames(all.phynotype) <-all.phynotype$sample  # Set sample as row names
all.phynotype <- all.phynotype[, -which(names(all.phynotype) == "sample")]  # Remove 'sample' column

################################################################################
# Align datasets by matching common sample IDs
sample_ids_gene <- colnames(all.gene_filtered)
sample_ids_meth <- colnames(all.meth_filtered)

# Find common sample IDs
common_samples <- intersect(colnames(all.gene_filtered), colnames(all.meth_filtered))

# Print the number of common samples
length(common_samples)

# Subset both datasets using the common samples
all.gene_filtered <- all.gene_filtered[, common_samples, drop = FALSE]
all.meth_filtered <- all.meth_filtered[, common_samples, drop = FALSE]

# Check the number of rows after filtering
nrow(all.gene_filtered)  # Should be > 0
nrow(all.meth_filtered)  # Should be > 0

########################################################
# Map Ensembl IDs to gene symbols and handle duplicates
ensembl_ids_no_version <- sub("\\..*", "", rownames(all.gene_filtered))  # Remove version suffix from Ensembl IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids_no_version, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Identify and handle duplicated gene symbols
duplicated_genes <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
gene_symbols[duplicated_genes] <- paste0(gene_symbols[duplicated_genes], "_", which(duplicated_genes))

# Ensure gene symbols are unique
length(unique(gene_symbols)) == length(gene_symbols)  # Should return TRUE

# Assign the unique gene symbols as row names
rownames(all.gene_filtered) <- gene_symbols

#################################################################
# Prepare phenotype groups for analysis
groups <- all.phynotype$sample_type  # Extract sample type as groups for analysis
t_all.meth_filtered=t(all.meth_filtered)  # Transpose methylation data
t_all.gene_filtered=t(all.gene_filtered)  # Transpose gene expression data

################################################################
# Perform PLS-DA with cross-validation
X <- cbind(t_all.meth_filtered, t_all.gene_filtered)  # Combine methylation and gene data
Y <- as.factor(groups)  # Convert groups to a factor for classification

# Fit PLS-DA model with 10 components
plsda_model <- plsda(X, Y, ncomp = 10)  

# Performance evaluation with cross-validation
set.seed(30)  # For reproducibility
perf_plsda <- perf(plsda_model, validation = 'Mfold', folds = 5, 
                   progressBar = FALSE, nrepeat = 10)

# Plot performance results
plot(perf_plsda, sd = TRUE, legend.position = 'horizontal')

###############################################################################
# Perform PCA on gene expression data
pca.gene <- pca(t_all.gene_filtered, ncomp = 4)
plotIndiv(pca.gene, group = groups, legend = TRUE, title = "Gene Expression PCA")

# Perform PCA on methylation data
pca.meth <- pca(t_all.meth_filtered, ncomp = 4)
plotIndiv(pca.meth, group = groups, legend = TRUE, title = "Methylation Expression PCA")

# Perform PLS integration between gene and methylation data
pls.result_1 <- pls(X = t_all.gene_filtered, Y = t_all.meth_filtered, ncomp = 4)
plotIndiv(pls.result_1, group = groups, legend = TRUE, title = "PLS Integration")

# Perform sparse PLS (sPLS) analysis
spls.result <- spls(X = t_all.gene_filtered, Y = t_all.gene_filtered, ncomp = 4, keepX = c(50, 50), keepY = c(50, 50))
plotVar(spls.result, cutoff = 0.8, title = "sPLS Correlation Circle", cex = c(6, 6))

# Perform PLS-DA on gene expression data
plsda.result <- plsda(X = t_all.gene_filtered, Y = groups, ncomp = 4)
plotIndiv(plsda.result, group = groups, legend = TRUE, title = "PLS-DA on Gene Expression")

# Perform DIABLO (Data Integration Analysis for Biomarker discovery using Latent cOmponents)
diablo.result <- block.plsda(X = list(gene = t_all.gene_filtered, meth = t_all.meth_filtered), Y = groups, ncomp = 4)
plotIndiv(diablo.result, legend = TRUE, title = "DIABLO Integration")

# Plot DIABLO variable plot
plotVar(diablo.result, cutoff = 0.85, title = "DIABLO Variable Plot", cex = c(6, 6))