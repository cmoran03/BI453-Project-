---
title: "BI453 Project"
author: "Caoimhe Moran"
date: "2026-03-07"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

This project applies Non-Negative Matrix Factorisation to miRNA samples of leukemia patients to cluster the different leukemia subtypes.


## Load Packages

```{r}
library(dplyr)
library(NMF)
library(mclust)
library(matrixcalc)
library(pheatmap)
library(ggplot2)
library(GEOquery)
```

## Pre-Processing

```{r}
# Remove rows and columns not needed for analysis
gset <- getGEO("GSE51908", GSEMatrix = TRUE)

#Create a matrix (miRNA x sample) with the expression data from the GEO dataset
GSE51908_series_matrix <- exprs(gset[[1]])

#Create a vector for the phenotypic data of the GEO dataset
pheno <- pData(gset[[1]])
head(pheno$title)

#Remove the control and cell line samples
A <- GSE51908_series_matrix[ ,-c(1:24, 43:78, 108:160, 174:190)]

# Make matrix non-negative
A <- nneg(A, method = "absolute", threshold = 0, shift = TRUE)
```

## Create Ground Truths

```{r}
head(colnames(A))

# Ground Truths when k=3
disease3 <- c("AML", "AML", "AML", "AML", "AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML",
"AML","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL")
covariates <- data.frame(Disease = disease3)
rownames(covariates) <- colnames(A)

# Ground Truths when k=2
disease2 <- c("AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL")
covariates2 <- data.frame(Disease = disease2)
rownames(covariates2) <- colnames(A)

# Ground Truths for ALL data
covariates_ALL <- data.frame(covariates[-c(1:18), ])
rownames(covariates_ALL) <- colnames(A[, -c(1:18)])
colnames(covariates_ALL) <- ("Disease")
```

## Rank Selection

```{r}
# Compare NMF for ranks k=2, k=3 and k=4 using the Brunet method
kcompar <- nmf(A, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

# Create consensus maps comparing each rank
png("RankComparison_kl.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar, annCol = covariates, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

# Create lines graphs comparing Cophenetic Correlation coefficient and silhouette of each rank
png("CopheneticCorrelation_Silhouette_RSS.png", width = 3000, height = 2400, res = 300)
plot(kcompar, what = c("cophenetic","rss"))
dev.off()

# Run NMF for each rank seperately
res_k2 <- nmf(A, rank = 2, method = "brunet", seed = "random", nrun = 50)
res_k3 <- nmf(A, rank = 3, method = "brunet", seed = "random", nrun = 50)
res_k4 <- nmf(A, rank = 4, method = "brunet", seed = "random", nrun = 50)

# Create cluster assignments based on NMF
clusters_k2<- predict(res_k2)
clusters_k3 <- predict(res_k3)
clusters_k4 <- predict(res_k4)

# Calculate ARI by comparing NMF clusters to Ground Truths to evaluate
# clustering accuracy
ari_k2 <- adjustedRandIndex(clusters_k2, covariates2$Disease)
ari_k3 <- adjustedRandIndex(clusters_k3, covariates$Disease)

# Run NMF for k=3 using the Lee method
res_k3_lee <- nmf(A, rank = 3, method = "lee", seed = "random", nrun = 50)

# Extract W and H matrices from NMF 
w_k3_lee <- basis(res_k3_lee)
h_k3_lee <- coef(res_k3_lee)

# Calculate residuals of NMF by subtracting the product of the W and H matrices from matrix A
residuals <- A - w_k3_lee%*%h_k3_lee

# Plot Residuals Quantiles against Theorotical Quantiles to visualise the distribution of the data
qqnorm(as.vector(residuals),
       main = "QQ Plot of NMF Residuals",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")

# Add line of best fit to the QQ Plot
qqline(as.vector(residuals), col = "red", lwd = 2)

```

## NMF by subtype

```{r}
# NMF for ALL data

# Assign ALL data as the 19th to 60th samples of A
ALL_data <- as.matrix(A[ ,(19:60)])

# Compare NMF for ranks k=2, k=3 and k=4 for ALL data
kcompar_ALL <- nmf(ALL_data, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

# Create consensus maps comparing each rank
png("RankSelection_kl_ALL.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar_ALL, annCol =  covariates_ALL, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

# Run NMF for best rank
resALL_k2 <- nmf(ALL_data, rank = 2, method = "brunet", seed = "random", nrun = 50)

# Create cluster assignments based on NMF
clustersALL_k2<- predict(resALL_k2)

# Calculate ARI by comparing NMF clusters to Ground Truths to evaluate clustering accuracy
ari_ALL_k2_kl <- adjustedRandIndex(clustersALL_k2, covariates_ALL$Disease)

# NMF for AML data

# Assign AML data as the 1st to 18th samples
AML_data <- as.matrix(A[ ,(1:18)])

# Compare NMF for ranks k=2, k=3 and k=4 for AML data
kcompar_AML <- nmf(AML_data, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

# Create consensus maps comparing each rank
png("RankSelection_kl_AML.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar_AML, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

# Create lines graphs comparing Cophenetic Correlation coefficient and silhouette of each rank
png("CopheneticCorrelation_Silhouette_RSS_AML.png", width = 3000, height = 2400, res = 300)
plot(kcompar_AML, what = c("cophenetic", "rss"))
dev.off()

# Run NMF for each rank seperately
resAML_k2 <- nmf(AML_data, rank = 2, method = "brunet", seed = "random", nrun = 50)
resAML_k3 <- nmf(AML_data, rank = 3, method = "brunet", seed = "random", nrun = 50)
resAML_k4 <- nmf(AML_data, rank = 4, method = "brunet", seed = "random", nrun = 50)

# Create cluster assignments based on NMF
clustersAML_k2 <- predict(resAML_k2)
clustersAML_k3 <- predict(resAML_k3)
clustersAML_k4 <- predict(resAML_k4)

# Calculate silhouette score to measure cluster seperation
sil_AML_k2 <- mean(silhouette(as.integer(clustersAML_k2), dist(t(AML_data)))[, 3])
sil_AML_k3 <- mean(silhouette(as.integer(clustersAML_k3), dist(t(AML_data)))[, 3])
sil_AML_k4 <- mean(silhouette(as.integer(clustersAML_k4), dist(t(AML_data)))[, 3])


```

## Post-Processing Max Normalisation

```{r}
# Extract W and H matrices from AML clustering
w_AML_k2 <- basis(resAML_k2)
h_AML_k2 <- coef(resAML_k2)

# Identify maximum value in each column of the W matrix 
scales_AML_k2 <- apply(w_AML_k2, 2, max)

# Divide each entry of each column in W by the maximum value in the column to max normalise W
w_AML_k2_norm <- sweep(w_AML_k2, 2, scales_AML_k2, "/")

# Multiply each entry of each row in H by the maimum value in the correspond column in W to max normalise H
h_AML_k2_norm <- sweep(h_AML_k2, 1, scales_AML_k2, "*")

# Assign each sample to the cluster corresponding to the row with the highest value  in each column of the max-normalised H matrix 
clusters_AML_k2_norm <- apply(h_AML_k2_norm, 2, which.max)

# Calculate absolute difference between cluster weights to identify miRNAs contributing most strongly to cluster separation
diff_w <- abs(w_AML_k2_norm[, 1] - w_AML_k2_norm[, 2])

# Select miRNAs with largest weight difference with a threshold that reduces dimensionality
top_miRNA <- AML_data[diff_w > 0.3, ]

# Run NMF for the cleaned AML data
AML_k2_cl <- nmf(top_miRNA, rank = 2, method = "brunet", seed = "random", nrun = 50)
AML_k3_cl <- nmf(top_miRNA, rank = 3, method = "brunet", seed = "random", nrun = 50)
AML_k4_cl <- nmf(top_miRNA, rank = 4, method = "brunet", seed = "random", nrun = 50)

# Create cluster assignments for cleaned AML data
clustersAML_k2_cl <- predict(AML_k2_cl)
clustersAML_k3_cl <- predict(AML_k3_cl)
clustersAML_k4_cl <- predict(AML_k4_cl)

# Calulate silhouette score for cleaned AML data
sil_AML_k2_cl <- mean(silhouette(as.integer(clustersAML_k2_cl), dist(t(top_miRNA)))[, 3])
sil_AML_k3_cl <- mean(silhouette(as.integer(clustersAML_k3_cl), dist(t(top_miRNA)))[, 3])
sil_AML_k4_cl <- mean(silhouette(as.integer(clustersAML_k4_cl), dist(t(top_miRNA)))[, 3])

# Store silhouette score of the AML data before and after filtering in a data frame
sil_df <- data.frame(
  Dataset = c(
    rep("AML", 3),
    rep("AML cleaned", 3)),
  k = rep(c(2, 3, 4), 2),
  Silhouette = c(
    sil_AML_k2, sil_AML_k3,
    sil_AML_k4, sil_AML_k2_cl,
    sil_AML_k3_cl, sil_AML_k4_cl))

# Visualise the difference of filtering on the silhouette score of the AML data
png("Silhouette_Normalisation.png")
ggplot(sil_df, aes(x = factor(k), y = Silhouette, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Dataset) +
  theme_minimal() +
  labs(
    title = "Silhouette Scores Before and After Removal of Low Variance miRNAs",
    x = "Number of Components (k)",
    y = "Average Silhouette Score" )
dev.off()

```

## Identitfy Differential miRNAs in AML clusters

```{r}
# Select top 20 miRNAs with the largest contribution differences
top20_wk2 <- sort(diff_w, decreasing = TRUE)[1:20]

# Store results in data frame for visualisation
top20_miRNAs_wk2 <- data.frame(miRNA = names(top20_wk2), Difference = top20_wk2)

# Visualise the top 20 miRNAs and their difference in contribution to the two clusters
png("Top_miRNAs_k2.png")
ggplot(top20_miRNAs_wk2, aes(
  x = reorder(miRNA, Difference), y = Difference, fill = Difference)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "blue") +
    labs(
     title = "Top 20 miRNAs Contributing to AML Cluster Separation when k=2",
     x = "miRNA",
     y = "Weight Difference")
dev.off()

```
