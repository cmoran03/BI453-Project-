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
#Remove rows and columns not needed for analysis
gset <- getGEO("GSE51908", GSEMatrix = TRUE)


GSE51908_series_matrix <- exprs(gset[[1]])

pheno <- pData(gset[[1]])
head(pheno$title)

A <- GSE51908_series_matrix[ ,-c(1:24, 43:78, 108:160, 174:190)]

#Make matrix non-negative
A <- nneg(A, method = "absolute", threshold = 0, shift = TRUE)
```

##Create Ground Truths

```{r}
head(colnames(A))

#Ground Truths when k=3
disease3 <- c("AML", "AML", "AML", "AML", "AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML",
"AML","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","B-cell ALL","T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL", "T-cell ALL")
covariates <- data.frame(Disease = disease3)
rownames(covariates) <- colnames(A)

#Ground Truths when k=2
disease2 <- c("AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","AML","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL",
"ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL")
covariates2 <- data.frame(Disease = disease2)
rownames(covariates2) <- colnames(A)
```

##Rank Selection

```{r}
kcompar <- nmf(A, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

png("RankComparison_kl.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar, annCol = covariates, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

png("CopheneticCorrelation_Silhouette_RSS.png", width = 3000, height = 2400, res = 300)
plot(kcompar, what = c("cophenetic","rss"))
dev.off()

res_k2 <- nmf(A, rank = 2, method = "brunet", seed = "random", nrun = 50)

clusters_k2<- predict(res_k2)
True_labels2 <- covariates2$Disease

ari_k2 <- adjustedRandIndex(clusters_k2, True_labels2)
sil_k2 <- mean(silhouette(as.integer(clusters_k2), dist(t(A)))[, 3])

res_k3 <- nmf(A, rank = 3, method = "brunet", seed = "random", nrun = 50)

clusters_k3 <- predict(res_k3)
True_labels3 <- covariates$Disease

ari_k3 <- adjustedRandIndex(clusters_k3, True_labels3)
sil_k3 <- mean(silhouette(as.integer(clusters_k3), dist(t(A)))[, 3])

res_k4 <- nmf(A, rank = 4, method = "brunet", seed = "random", nrun = 50)
clusters_k4 <- predict(res_k4)
sil_k4 <- mean(silhouette(as.integer(clusters_k4), dist(t(A)))[, 3])

sil_comp <- data.frame(Silhouette = c(sil_k2, sil_k3, sil_k4), 
                       Rank = c(2, 3, 4))

png("Silhouette_Rank_Compare.png")
ggplot(sil_comp, aes(x = factor(Rank), y = Silhouette, group = 1)) +
  geom_line(stat = "identity") +
  geom_point(size = 3) +
  labs(
    x = "Rank (k)",
    y = "Silhouette Score",
    title = "Silhouette Scores Across NMF Ranks")
dev.off()


kcompar_lee <- nmf(A, rank = 2:4, method = "lee", seed = "random", nrun = 50)

png("RankComparison_lee.pdf", width = 3000, height = 2400, res = 300)
consensusmap(kcompar_lee, annCol = covariates, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

png("CopheneticCorrelation_lee.png", width = 3000, height = 2400, res = 300)
plot(kcompar_lee, what = c("cophenetic", "rss"))
dev.off()


res_k3_lee <- nmf(A, rank = 3, method = "lee", seed = "random", nrun = 50)
w_k3_lee <- basis(res_k3_lee)
h_k3_lee <- coef(res_k3_lee)

residuals <- A - w_k3_lee%*%h_k3_lee
qqnorm(as.vector(residuals),
       main = "QQ Plot of NMF Residuals",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(as.vector(residuals), col = "red", lwd = 2)

```

##NMF by subtype

```{r}

ALL_data <- as.matrix(A[ ,(19:60)])

kcompar_ALL <- nmf(ALL_data, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

covariates_ALL <- data.frame(covariates[-c(1:18), ])
rownames(covariates_ALL) <- colnames(A[, -c(1:18)])
colnames(covariates_ALL) <- ("Disease")

png("RankSelection_kl_ALL.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar_ALL, annCol =  covariates_ALL, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

resALL_k2 <- nmf(ALL_data, rank = 2, method = "brunet", seed = "random", nrun = 50)
clustersALL_k2<- predict(resALL_k2)
True_labels_ALL <- covariates_ALL$Disease

ari_ALL_k2_kl <- adjustedRandIndex(clustersALL_k2, True_labels_ALL)


AML_data <- as.matrix(A[ ,(1:18)])

kcompar_AML <- nmf(AML_data, rank = 2:4, method = "brunet", seed = "random", nrun = 50)

png("RankSelection_kl_AML.png", width = 3000, height = 2400, res = 300)
consensusmap(kcompar_AML, labRow = NA, labCol = NA, tracks = c("silhouette"))
dev.off()

png("CopheneticCorrelation_Silhouette_RSS_AML.png", width = 3000, height = 2400, res = 300)
plot(kcompar_AML, what = c("cophenetic", "rss"))
dev.off()

resAML_k2 <- nmf(AML_data, rank = 2, method = "brunet", seed = "random", nrun = 50)
clustersAML_k2 <- predict(resAML_k2)
sil_AML_k2 <- mean(silhouette(as.integer(clustersAML_k2), dist(t(AML_data)))[, 3])

resAML_k3 <- nmf(AML_data, rank = 3, method = "brunet", seed = "random", nrun = 50)
clustersAML_k3 <- predict(resAML_k3)
sil_AML_k3 <- mean(silhouette(as.integer(clustersAML_k3), dist(t(AML_data)))[, 3])

resAML_k4 <- nmf(AML_data, rank = 4, method = "brunet", seed = "random", nrun = 50)
clustersAML_k4 <- predict(resAML_k4)
sil_AML_k4 <- mean(silhouette(as.integer(clustersAML_k4), dist(t(AML_data)))[, 3])


```

## Post-Processing Max Normalisation

```{r}
w_AML_k2 <- basis(resAML_k2)
h_AML_k2 <- coef(resAML_k2)

scales_AML_k2 <- apply(w_AML_k2, 2, max)
w_AML_k2_norm <- sweep(w_AML_k2, 2, scales_AML_k2, "/")
h_AML_k2_norm <- sweep(h_AML_k2, 1, scales_AML_k2, "*")
clusters_AML_k2_norm <- apply(h_AML_k2_norm, 2, which.max)

diff_w <- abs(w_AML_k2_norm[, 1] - w_AML_k2_norm[, 2])

top_miRNA <- AML_data[diff_w > 0.3, ]

AML_k2_cl <- nmf(top_miRNA, rank = 2, method = "brunet", seed = "random", nrun = 50)
clustersAML_k2_cl <- predict(AML_k2_cl)

AML_k3_cl <- nmf(top_miRNA, rank = 3, method = "brunet", seed = "random", nrun = 50)
clustersAML_k3_cl <- predict(AML_k3_cl)

AML_k4_cl <- nmf(top_miRNA, rank = 4, method = "brunet", seed = "random", nrun = 50)
clustersAML_k4_cl <- predict(AML_k4_cl)

sil_AML_k2_cl <- mean(silhouette(as.integer(clustersAML_k2_cl), dist(t(top_miRNA)))[, 3])
sil_AML_k3_cl <- mean(silhouette(as.integer(clustersAML_k3_cl), dist(t(top_miRNA)))[, 3])
sil_AML_k4_cl <- mean(silhouette(as.integer(clustersAML_k4_cl), dist(t(top_miRNA)))[, 3])


sil_df <- data.frame(
  Dataset = c(
    rep("AML", 3),
    rep("AML cleaned", 3)),
  k = rep(c(2, 3, 4), 2),
  Silhouette = c(
    sil_AML_k2, sil_AML_k3,
    sil_AML_k4, sil_AML_k2_cl,
    sil_AML_k3_cl, sil_AML_k4_cl))

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

##Identitfy Differential miRNAs in AML clusters

```{r}
top20_wk2 <- sort(diff_w, decreasing = TRUE)[1:20]
top20_miRNAs_wk2 <- data.frame(miRNA = names(top20_wk2), Difference = top20_wk2)

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





Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
