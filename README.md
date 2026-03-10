title: "BI453 Project"
author: "Caoimhe Moran"
date: "2026-03-07"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction

This project applies Non-Negative MAtrix Factorisation to miRNA samples of leukemia patients to cluster the different leukemia subtypes.


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
Leuk_data <- GSE51908_series_matrix[-c(1:27, 29:63, 911), ]
Leuk_data <- Leuk_data[ ,-c(2:25, 44:79, 109:161, 175:191)]

#Convert to data frame
Leuk_data <- as.data.frame(Leuk_data)
          
#Convert row names to miRNA names                  
rownames(Leuk_data) <- Leuk_data[, 1]
Leuk_data <- Leuk_data[, -1]

#Convert column names to sample names
colnames(Leuk_data) <- Leuk_data[ 1,]
Leuk_data <- Leuk_data[ -1,]

#Convert data frame to matrix
A <- as.matrix(sapply(Leuk_data, as.numeric))
rownames(A) <- rownames(Leuk_data)

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

##NMF 

```{r}
res_k3 <- nmf(A, rank = 3, method = "lee", seed = "random", nrun = 50)
w_k3 <- basis(res_k3)
dim(w_k3)
h_k3 <- coef(res_k3)
dim(h_k3)

residuals <- A - w_k3%*%h_k3
qqnorm(as.vector(residuals))
qqline(as.vector(residuals))

clusters_k3 <- predict(res_k3)
True_labels3 <- covariates$Disease

ari_k3 <- adjustedRandIndex(clusters_k3, True_labels3)

sil_k3 <- mean(silhouette(as.integer(clusters_k3), dist(t(A)))[, 3])

res_k2 <- nmf(A, rank = 2, method = "lee", seed = "random", nrun = 50)
w_k2 <- basis(res_k2)
h_k2 <- coef(res_k2)
consensusmap(res_k2)

clusters_k2<- predict(res_k2)
True_labels2 <- covariates2$Disease

ari_k2 <- adjustedRandIndex(clusters_k2, True_labels2)

sil_k2 <- mean(silhouette(as.integer(clusters_k2), dist(t(A)))[, 3])

res_k4 <- nmf(A, rank = 4, method = "lee", seed = "random", nrun = 50)
clusters_k4 <- predict(res_k4)
sil_k4 <- mean(silhouette(as.integer(clusters_k4), dist(t(A)))[, 3])

```

##Rank Selection

```{r}

kcompar <- nmf(A, rank = 2:4, method = "lee", seed = "random", nrun = 10)

pdf("RankComparison_lee.pdf", width = 14, height = 12)
consensusmap(kcompar, annCol = covariates)
dev.off()

pdf("CopheneticCorrelation_lee.pdf", width = 14, height = 12)
plot(kcompar, what = c("cophenetic"))
dev.off()

sil_comp <- data.frame(Silhouette = c(sil_k2, sil_k3, sil_k4), 
                       Rank = c(2, 3, 4))

pdf("Silhouette_Rank_Compare.pdf")
ggplot(sil_comp, aes(x = factor(Rank), y = Silhouette, group = 1)) +
  geom_line(stat = "identity") +
  geom_point(size = 3) +
  labs(
    x = "Rank (k)",
    y = "Silhouette Score",
    title = "Silhouette Scores Across NMF Ranks")
dev.off()

```

##NMF by subtype

```{r}

AML_data <- as.matrix(A[ ,(1:18)])

ALL_data <- as.matrix(A[ ,(19:60)])

resAML_k2 <- nmf(AML_data, rank = 2, method = "lee", seed = "random", nrun = 50)
w_AML_k2 <- basis(resAML_k2)
row.names(w_AML_k2) <- row.names(A)
h_AML_k2 <- coef(resAML_k2)
pdf("AML_k2_lee.pdf")
consensusmap(resAML_k2)
dev.off()

clustersAML_k2 <- predict(resAML_k2)
sil_AML_k2 <- mean(silhouette(as.integer(clustersAML_k2), dist(t(AML_data)))[, 3])

resAML_k3 <- nmf(AML_data, rank = 3, method = "lee", seed = "random", nrun = 50)
w_AML_k3 <- basis(resAML_k3)
h_AML_k3 <- coef(resAML_k3)
pdf("AML_k3_lee.pdf")
consensusmap(resAML_k3)
dev.off()

clustersAML_k3 <- predict(resAML_k3)
sil_AML_k3 <- mean(silhouette(as.integer(clustersAML_k3), dist(t(AML_data)))[, 3])

covariates_ALL <- data.frame(covariates[-c(1:18), ])
rownames(covariates_ALL) <- colnames(A[, -c(1:18)])
colnames(covariates_ALL) <- ("Disease")

resALL_k2 <- nmf(ALL_data, rank = 2, method = "lee", seed = "random", nrun = 50)
w_ALL_k2 <- basis(resALL_k2)
h_ALL_k2 <- coef(resALL_k2)
pdf("ALL_k2_lee.pdf")
consensusmap(resALL_k2, annCol = covariates_ALL)
dev.off()

clustersALL_k2<- predict(resALL_k2)
table(clustersALL_k2)
head(clustersALL_k2)
True_labels_ALL <- covariates_ALL$Disease
table(True_labels_ALL)
head(True_labels_ALL)

ari_ALL_k2 <- adjustedRandIndex(clustersALL_k2, True_labels_ALL)

sil_ALL_k2 <- mean(silhouette(as.integer(clustersALL_k2), dist(t(ALL_data)))[, 3])

resALL_k3 <- nmf(ALL_data, rank = 3, method = "lee", seed = "random", nrun = 50)
w_ALL_k3 <- basis(resALL_k3)
h_ALL_k3 <- coef(resALL_k3)
pdf("ALL_k3_lee.pdf")
consensusmap(resALL_k3)
dev.off()

clustersALL_k3<- predict(resALL_k3)
table(clustersALL_k3)
head(clustersALL_k3)

ari_ALL_k3 <- adjustedRandIndex(clustersALL_k3, True_labels_ALL)

sil_ALL_k3 <- mean(silhouette(as.integer(clustersALL_k3), dist(t(ALL_data)))[, 3])

```

## Post-Processing Normalisation

```{r}
scales_k3 <- apply(w_k3, 2, max)
w_k3_norm <- sweep(w_k3, 2, scales_k3, "/")
h_k3_norm <- sweep(h_k3, 1, scales_k3, "*")
clusters_k3_norm <- apply(h_k3_norm, 2, which.max)
ari_k3_norm <- adjustedRandIndex(clusters_k3_norm, True_labels3)

scales_k2 <- apply(w_k2, 2, max)
w_k2_norm <- sweep(w_k2, 2, scales_k2, "/")
h_k2_norm <- sweep(h_k2, 1, scales_k2, "*")
clusters_k2_norm <- apply(h_k2_norm, 2, which.max)
ari_k2_norm <- adjustedRandIndex(clusters_k2_norm, True_labels2)

scales_ALL_k2 <- apply(w_ALL_k2, 2, max)
w_ALL_k2_norm <- sweep(w_ALL_k2, 2, scales_ALL_k2, "/")
h_ALL_k2_norm <- sweep(h_ALL_k2, 1, scales_ALL_k2, "*")
clusters_ALL_k2_norm <- apply(h_ALL_k2_norm, 2, which.max)
ari_ALL_k2_norm <- adjustedRandIndex(clusters_ALL_k2_norm, True_labels_ALL)

scales_ALL_k3 <- apply(w_ALL_k3, 2, max)
w_ALL_k3_norm <- sweep(w_ALL_k3, 2, scales_ALL_k3, "/")
h_ALL_k3_norm <- sweep(h_ALL_k3, 1, scales_ALL_k3, "*")
clusters_ALL_k3_norm <- apply(h_ALL_k3_norm, 2, which.max)
ari_ALL_k3_norm <- adjustedRandIndex(clusters_ALL_k3_norm, True_labels_ALL)

scales_AML_k2 <- apply(w_AML_k2, 2, max)
w_AML_k2_norm <- sweep(w_AML_k2, 2, scales_AML_k2, "/")
h_AML_k2_norm <- sweep(h_AML_k2, 1, scales_AML_k2, "*")
clusters_AML_k2_norm <- apply(h_AML_k2_norm, 2, which.max)

scales_AML_k3 <- apply(w_AML_k3, 2, max)
w_AML_k3_norm <- sweep(w_AML_k3, 2, scales_AML_k3, "/")
h_AML_k3_norm <- sweep(h_AML_k3, 1, scales_AML_k3, "*")
clusters_AML_k3_norm <- apply(h_AML_k3_norm, 2, which.max)

sil_AML_k3_norm <- mean(silhouette(clusters_AML_k3_norm, dist(t(AML_data)))[, 3])
sil_AML_k2_norm <- mean(silhouette(clusters_AML_k2_norm, dist(t(AML_data)))[, 3])
sil_ALL_k3_norm <- mean(silhouette(clusters_ALL_k3_norm, dist(t(ALL_data)))[, 3])
sil_ALL_k2_norm <- mean(silhouette(clusters_ALL_k2_norm, dist(t(ALL_data)))[, 3])
sil_k2_norm <- mean(silhouette(clusters_k2_norm, dist(t(A)))[, 3])
sil_k3_norm <- mean(silhouette(clusters_k3_norm, dist(t(A)))[, 3])

sil_df <- data.frame(
  Dataset = c(
    rep("ALL", 4),
    rep("AML", 4),
    rep("Full", 4)),
  k = rep(c(2, 2, 3, 3), 3),
  Normalisation = rep(c("Pre", "Post", "Pre", "Post"), 3),
  Silhouette = c(
    sil_ALL_k2, sil_ALL_k2_norm,
    sil_ALL_k3, sil_ALL_k3_norm,
    sil_AML_k2, sil_AML_k2_norm,
    sil_AML_k3, sil_AML_k3_norm,
    sil_k2, sil_k2_norm,
    sil_k3, sil_k3_norm ))

pdf("Silhouette_Normalisation.pdf")
ggplot(sil_df, aes(x = factor(k), y = Silhouette, fill = Normalisation)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Dataset) +
  theme_minimal() +
  labs(
    title = "Silhouette Scores Before and After Max Normalisation",
    x = "Number of Components (k)",
    y = "Average Silhouette Score" )
dev.off()

ari_df <- data.frame(
  Dataset = c (
    rep("ALL",4), 
    rep("Full", 4)),
  k = rep(c(2,2,3,3), 2),
  Normalisation = rep(c("Pre", "Post", "Pre", "Post"), 2),
  ARI = c(
    ari_ALL_k2, ari_ALL_k2_norm,
    ari_ALL_k3, ari_ALL_k3_norm,
    ari_k2, ari_k2_norm,
    ari_k3, ari_k3_norm)) 
    
pdf("ARI_Normalisation.pdf")
ggplot(ari_df, aes(x = factor(k), y = ARI, fill = Normalisation)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Dataset) +
  theme_minimal() +
  labs(
    title = "ARI Scores Before and After Max Normalisation",
    x = "Number of Components (k)",
    y = "Average ARI Score" )
dev.off()

```

##Identitfy Differential miRNAs in AML clusters

```{r}
diff_w <- abs(w_AML_k2_norm[, 1] - w_AML_k2_norm[, 2])
top20_wk2 <- sort(diff_w, decreasing = TRUE)[1:20]
top20_miRNAs_wk2 <- data.frame(miRNA = names(top20_wk2), Difference = top20_wk2)

pdf("Top_miRNAs_k2.pdf")
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

diff_w3 <- apply(w_AML_k3_norm, 1, function(x) max(x) - min(x))
top20_wk3 <- sort(diff_w3, decreasing = TRUE)[1:20]
top20_miRNAs_wk3 <- data.frame(miRNA = names(top20_wk3), Difference = top20_wk3)

pdf("Top_miRNAs_k3.pdf")
ggplot(top20_miRNAs_wk3, aes(
  x = reorder(miRNA, Difference), y = Difference, fill = Difference)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "blue") +
  labs(
    title = "Top 20 miRNAs Contributing to AML Cluster Separation when k=3",
    x = "miRNA",
    y = "Weight Difference")
dev.off()





Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
