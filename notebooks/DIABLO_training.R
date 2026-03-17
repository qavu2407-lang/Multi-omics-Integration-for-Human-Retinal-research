# ==================== Install and Load Packages ====================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("mixOmics", quietly = TRUE))
  BiocManager::install("mixOmics")

library(mixOmics)
library(ggplot2)
library(dplyr)
set.seed(123) # for reproducibility

# ==================== Load Data ====================
cat("Loading data...\n")

# Import data and inspect it
head(proteins)
head(lipids)
head(metabolites)

