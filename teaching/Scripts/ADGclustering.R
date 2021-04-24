
### Required packages

## Installation of FactoMineR (only if needed)

install.packages("FactoMineR")
require(FactoMineR)

## Installation of mvtnorm and fdrtool (only if needed)

install.packages("flashClust")
install.packages("gplots")

require(flashClust)
require(gplots)

### Set working directory

setwd("C:/Users/David/Dropbox/ADG2020/Slides/Session4/Data")
dir()                                                

### Import data

## Import gene expression data
expressions = read.table("foiePositive.txt",header=TRUE)
dim(expressions)           # Numbers of rows and columns
head(expressions[,1:10])   # Displays first 6 rows of first 10 columns 

## Import heterogeneity-adjusted gene expression data
adjusted_expressions = read.table("foiePositiveAdjusted.txt",header=TRUE)
dim(adjusted_expressions)           # Numbers of rows and columns
head(adjusted_expressions[,1:10])   # Displays first 6 rows of first 10 columns 

## Import experimental design variables
covariates = read.table("covariates.txt",header=TRUE)
dim(covariates)            # Numbers of rows and columns
str(covariates)            # Overview of the data

## Import supplementary variables on chickens
external = read.table("external.txt",header=TRUE)
dim(external)            # Numbers of rows and columns
str(external)            # Overview of the data

### Clustering

# Clustering microarrays

x = as.matrix(expressions)
heatmap(x, Rowv = FALSE, col = rev(heat.colors(16))) 

hc = flashClust::hclust(dist(t(x)),method="ward")
plot(hc,frame.plot=TRUE,labels=NULL,xlab="Clusters",main="Ward's cluster dendogram",
     sub="Distance: Euclidean",
     cex.axis=1.25,cex.main=1.25,cex.lab=1.25)
hc.cut = cutree(hc,k=3)

# Spectral clustering (clustering on principal components)

# Create a data table with expression data and supplmentary variables
dta = data.frame(covariates,PectMajor_g=external$PectMajor_g,t(x))
head(dta)

# PCA of expresssion data with summplementary variables
dta.PCA = PCA(dta,quali.sup=1:2,quanti.sup=3)
# Eigenvalues (useful to choose the number of PCs)
dta.PCA$eig               
# Correlations between supplementary covariate and PCs
dta.PCA$quanti.sup$cor
# Association between supplementary factors and PCs
dta.PCA$quali.sup$v.test  

# Spectral clustering with 10 PCs
dta.PCA = PCA(dta,quali.sup=1:2,quanti.sup=3,ncp=10)
foie.HCPC = HCPC(dta.PCA,consol=TRUE,proba=1e-06)

# Distribution in clusters
clusters = foie.HCPC$data.clust$clust
table(clusters,dta$Diet)

# Significant associations between clusters and supplementary factors
foie.HCPC$desc.var$test.chi2
foie.HCPC$desc.var$category

# Significant associations between clusters and supplementary covariate
head(foie.HCPC$desc.var$quanti.var)
# Significant associations between clusters and gene expressions (focus on cluster 1)
desc_cluster1 = foie.HCPC$desc.var$quanti[[1]]
head(desc_cluster1)
negative_associations_1 = rownames(desc_cluster1)[desc_cluster1[,1]<0] 
negative_associations_1 = is.element(rownames(x),negative_associations_1) 

# Clustering genes

## Heatmap of the correlation matrix
hc = heatmap.2(cor(t(x)), symm = TRUE, distfun = function(c) as.dist(sqrt(1 - c^2)), 
               trace="none",dendrogram="col")

## Ward's clustering of the correlation matrix 
hc = flashClust::hclust(as.dist(sqrt(1-cor(t(x))^2)),method="ward")
plot(hc,frame.plot=TRUE,labels=FALSE,xlab="Clusters",main="Ward's cluster dendogram",
     sub="Squared Distance: 1-squared correlation",
     cex.axis=1.25,cex.main=1.25,cex.lab=1.25)

# Cut the dendogram
gene_clusters = cutree(hc,k=2)
table(gene_clusters)

# Relationship with the results of microarray clustering
table(gene_clusters,negative_associations_1)

# Relationship with plasma composition 
external_plasma = external[,12:20]
cor_cluster2_external = cor(t(x)[,gene_clusters==2],external_plasma,
                            use="pairwise.complete.obs")

# Maximum absolute correlations 
colnames(external_plasma)[apply(abs(cor_cluster2_external),1,which.max)]

cor_cluster2_butyrate = cor(t(x)[,gene_clusters==2],
                      external_plasma[,c("plasma_B.OH")],
                      use="pairwise.complete.obs")
summary(cor_cluster2_butyrate)

plot(external$"plasma_B.OH",expressions["lip_FASN",],
     xlab="BOH",ylab="Expression of lip_FASN",
     main="lip_FASN expression",pch=16,bty="l",cex.axis=1.25,
     col=ifelse(covariates$Diet=="HL","orange","darkgray"),cex.lab=1.25)
legend("bottomright",bty="n",pch=16,cex=1.25,
       col=c("orange","darkgray"),legend=c("Diet: HL","Diet: BL"))
