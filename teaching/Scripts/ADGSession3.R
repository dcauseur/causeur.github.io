
### Required packages

## Installation of limma (only if needed)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma",update=FALSE)

require(limma)

## Installation of FAMT (only if needed)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute",update=FALSE)

install.packages("FAMT")

require(FAMT)

## Installation of fdrtool (only if needed)

install.packages("fdrtool")

require(fdrtool)

### Set working directory

setwd("C:/Users/David/Dropbox/ADG2020/Slides/Session3/Data")
dir()                                                

### Import data

## Import gene expression data
expressions = read.table("foie.txt",header=TRUE)
dim(expressions)           # Numbers of rows and columns
head(expressions[,1:10])   # Displays first 6 rows of first 10 columns 

## Import experimental design variables
covariates = read.table("covariates.txt",header=TRUE)
dim(covariates)            # Numbers of rows and columns
str(covariates)            # Overview of the data

### Choice of the test statistic

## Linear modeling using FAMT

# First, create an FAMTdata object gathering expressions and covariates
# Column names of expressions should correspond to an ID variable in covariates
# An ID variable is added in 1st column of covariates
covariates = data.frame(ID=rownames(covariates),covariates)

# as.FAMT data creates the FAMTdata object - idcovar gives the column number of the
# ID variable in covariates
foie.famt = as.FAMTdata(expressions,covariates,idcovar=1)

# Overview of data in the FAMTdata object
summaryFAMT(foie.famt)

# Fits the linear model: here, Diet+Genotype (x = 2nd and 3rd columns of covariates)
# Test for the Diet effect (test = 2nd column of covariates)
# nbf = 0 for a standard fit
foie.lmfit = modelFAMT(foie.famt,x=c(2,3),test=2,nbf=0)

# Display the names of the 9 components in foie.lmfit
names(foie.lmfit)

# p-values for the diet effect can be found in $pval
pvalD = foie.lmfit$pval
hist(pvalD,proba=TRUE,col="orange",xlab="p-values",
     main="Diet effect on gene expression",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)

# Genes with p-values <= 0.05
positives = pvalD<=0.05       # Identifies positive genes
mean(positives)               # Proportion of positive genes 
sort(pvalD[positives])[1:10]  # First 10 most significant genes

### Simulation of null and nonnull gene expressions

## A gene expression dataset with no signal 

n = ncol(expressions)                  # Number of microarrays                              
m = nrow(expressions)                  # Number of genes                         

# Generate a m x n random matrix of expression data
dta0 = matrix(rnorm(n*m),nrow=m,ncol=n)
# Give same names to columns of dta0 as in expressions
colnames(dta0) = colnames(expressions)

# Calculate p-values for group effect as above with FAMT

# as.FAMT data creates the FAMTdata object - idcovar gives the column number of the
# ID variable in covariates
dta0.famt = as.FAMTdata(dta0,covariates,idcovar=1)

# Fits the linear model: here, Diet+Genotype (x = 2nd and 3rd columns of covariates)
# Test for the Diet effect (test = 2nd column of covariates)
# nbf = 0 for a standard fit
dta0.lmfit = modelFAMT(dta0.famt,x=c(2,3),test=2,nbf=0)

# p-values for the diet effect can be found in $pval
pval0 = dta0.lmfit$pval
hist(pval0,proba=TRUE,nclass=20,col="orange",
     main="Histogram of p-values",xlab="p-values")
mtext("Data simulated with no diet effect")

# Genes with p-values <= 0.05
positives0 = pval0<=0.05        # Identifies positive genes
mean(positives0)                # Proportion of positive genes 
sort(pval0[positives0])[1:10]   # First 10 most significant genes

## A gene expression dataset with 5% of nonnull genes

# First, let us start from the dataset with no diet effect
dta1 = dta0    

# Then, let us set a diet effect on the first 2303 genes (5% of 46060)
# For these genes, the mean for diet=HL is 1, the mean for diet=BL remains 0 
dta1[1:2303,covariates$Diet=="HL"] = dta1[1:2303,covariates$Diet=="HL"] + 1 

# Calculate p-values for group effect as above with FAMT

# as.FAMT data creates the FAMTdata object - idcovar gives the column number of the
# ID variable in covariates
dta1.famt = as.FAMTdata(dta1,covariates,idcovar=1)

# Fits the linear model: here, Diet+Genotype (x = 2nd and 3rd columns of covariates)
# Test for the Diet effect (test = 2nd column of covariates)
# nbf = 0 for a standard fit
dta1.lmfit = modelFAMT(dta1.famt,x=c(2,3),test=2,nbf=0)

# p-values for the diet effect can be found in $pval
pval1 = dta1.lmfit$pval
hist(pval1,proba=TRUE,nclass=20,col="orange",
     main="Histogram of p-values",xlab="p-values")
mtext("Data simulated with diet effect on 5% of genes") 

# Genes with p-values <= 0.05
positives1 = pval1<=0.05        # Identifies positive genes
mean(positives1)                # Proportion of positive genes 
sort(pval1[positives1])[1:10]   # First 10 most significant genes

# Number of true positives
sum(pval1[1:2303]<=0.05) 

# Number of false positives
sum(pval1[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(pval1[-(1:2303)]<=0.05)/sum(pval1<=0.05)     

## Bonferroni correction

Bonfpval1 = p.adjust(pval1,method="bonferroni") # adjusted p-values (Bonferroni)
BonfPositives1 = Bonfpval1<=0.05 

# Number of true positives
sum(Bonfpval1[1:2303]<=0.05) 

# Number of false positives
sum(Bonfpval1[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(Bonfpval1[-(1:2303)]<=0.05)/sum(Bonfpval1<=0.05)     

## Benjamini-Hochberg correction

BHpval1 = p.adjust(pval1,method="BH") # adjusted p-values (BH)
BHPositives1 = BHpval1<=0.05 

# Number of true positives
sum(BHpval1[1:2303]<=0.05) 

# Number of false positives
sum(BHpval1[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(BHpval1[-(1:2303)]<=0.05)/sum(BHpval1<=0.05)     

## q-values

# On simulated data

pi0 = pval.estimate.eta0(pval1)   # Estimation of the proportion of true nulls
qvalues1 = pi0*BHpval1            # q-values

# Number of true positives
sum(qvalues1[1:2303]<=0.05) 

# Number of false positives
sum(qvalues1[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(qvalues1[-(1:2303)]<=0.05)/sum(qvalues1<=0.05)     

# On the expression dataset

pi0 = pval.estimate.eta0(pvalD)       # Estimation of the proportion of true nulls
BHpvalD = p.adjust(pvalD,method="BH") # adjusted p-values (BH)
qvaluesD = pi0*BHpvalD                # q-values

positives = qvaluesD<=0.05

# Number of positive genes with a control of the FDR at level 0.05
sum(positives)                         

# Suppose we choose to keep 100 genes: what is the estimated FDR?
sort(qvaluesD)[100]

## Moderated tests in limma

# First, set the design matrix
design = model.matrix(~Genotype+Diet,data=covariates)   
head(design)

# Then, fit the 2-way analysis of variance model
fit = lmFit(expressions,design)
head(fit$coefficients)          # Display the first 6 rows

# Now, calculate the moderated tests statistics and corresponding p-values
fit = eBayes(fit)      

# Display the top 10 most significant genes
topTable(fit, coef=3)  # Coef=column number in fit$coefficients for the test   

# Summarize the decisions
results = decideTests(fit,adjust.method="BH",p.value=0.05) 
print(summary(results)) 

# Volcano plot
BHmoderatedpval = p.adjust(fit$p.value[,"DietHL"],method="BH")
logFC = fit$coefficients[,3]

plot(logFC,-log10(BHmoderatedpval),xlab = "log-Fold Change",
     ylab = "-log10(pvalue)",main = "Volcano plot",cex = 0.6, pch = 19)
points(logFC[BHmoderatedpval <= 0.05],
       -log10(BHmoderatedpval)[BHmoderatedpval <= 0.05],
       cex = 0.6, pch = 19, col = "red")
abline(h = -log10(0.05),col = "blue",lty =2,lwd = 1.5)

# Gene selection

# Condition to be selected: BH adjusted p-value<=0.05 + |logFC|>=1
select = (BHmoderatedpval<=0.05) & (abs(logFC)>=1)

# Number of selected genes
sum(select)

### Heterogeneity in gene expressions

# A simulated heterogeneous gene expression dataset

z = rnorm(n)                                 
dta2 = sweep(dta1,MARGIN=2,STATS=z,FUN="+")                     

# Calculate p-values for group effect as above with FAMT

# as.FAMT data creates the FAMTdata object - idcovar gives the column number of the
# ID variable in covariates
dta2.famt = as.FAMTdata(dta2,covariates,idcovar=1)

# Fits the linear model: here, Diet+Genotype (x = 2nd and 3rd columns of covariates)
# Test for the Diet effect (test = 2nd column of covariates)
# nbf = 0 for a standard fit
dta2.lmfit = modelFAMT(dta2.famt,x=c(2,3),test=2,nbf=0)

# p-values for the diet effect can be found in $pval
pval2 = dta2.lmfit$pval
hist(pval2,proba=TRUE,nclass=20,col="orange",main="Histogram of p-values",
     xlab="p-values")
mtext("Heterogeneous data")

# BH correction
BHpval2 = p.adjust(pval2,method="BH")

# Genes with p-values <= 0.05
positives2 = BHpval2<=0.05        # Identifies positive genes
mean(positives2)                # Proportion of positive genes 

# Number of true positives
sum(BHpval2[1:2303]<=0.05) 

# Number of false positives
sum(BHpval2[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(BHpval2[-(1:2303)]<=0.05)/sum(BHpval2<=0.05)     

## Identification of heterogeneity

# Number of heterogeneity components

nbf = nbfactors(dta2.famt,x=c(2,3),test=2,diagnostic.plot=TRUE)

# Extraction of heterogeneity components

dta2.lmfit = modelFAMT(dta2.famt,x=c(2,3),test=2,nbf=1)

hist(dta2.lmfit$adjpval,proba=TRUE,nclass=20,col="orange",
     main="Histogram of heterogeneity-adjusted p-values",
     xlab="p-values")
mtext("Heterogeneous data")

# Comparison of selection performance with and without heterogeneity adjustment
results = summaryFAMT(dta2.lmfit,alpha=seq(0,0.05,0.01),pi0=NULL)
results

# q-values based on heterogeneity-adjusted p-values 
qvalues2 = results$pi0*p.adjust(dta2.lmfit$adjpval,method="BH")

# Genes with q-values <= 0.05
positives2 = qvalues2<=0.05        # Identifies positive genes
mean(positives2)                   # Proportion of positive genes 

# Number of true positives
sum(qvalues2[1:2303]<=0.05) 

# Number of false positives
sum(qvalues2[-(1:2303)]<=0.05)                     

# False Discovery Proportion
sum(qvalues2[-(1:2303)]<=0.05)/sum(qvalues2<=0.05)     

## Modeling heterogeneity in expression data

foie.lmfit = modelFAMT(foie.famt,x=c(2,3),test=2,nbf=NULL)
results = summaryFAMT(foie.lmfit,alpha=0.05,pi0=NULL)
results$pi0
results$nbreject

# Heterogeneity-adjusted p-values for the diet effect can be found in $adjpval
pvalD = foie.lmfit$adjpval
hist(pvalD,proba=TRUE,col="orange",xlab="p-values",
     main="Diet effect on gene expression",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
mtext("After heterogeneity adjustment")

# Genes with p-values <= 0.05
qvaluesD = results$pi0*p.adjust(foie.lmfit$adjpval,method="BH")
positives = qvaluesD<=0.05       # Identifies positive genes
mean(positives)               # Proportion of positive genes 
sort(pvalD[positives])[1:10]  # First 10 most significant genes

### Export tables

## Choose your selection method

# Either moderated tests + BH correction + |logFC|>=1
select1 = (BHmoderatedpval<=0.05) & (abs(logFC)>=1)

# Number of selected genes
sum(select1)

# Or heterogeneity adjustment + q-values + |logFC|>=1
select2 = (qvaluesD<=0.05) & (abs(logFC)>=1)

# Number of selected genes
sum(select2)

# Brief comparison
table(select1,select2)

## Then export the expression data for selected genes

foiePositive = expressions[select1,]
write.table(foiePositive,"foiePositive.txt")

## You can also export the heterogeneity-adjusted expression data for selected genes

expressions_adjusted = foie.lmfit$adjdata$expression
foiePositive_adjusted = expressions_adjusted[select1,]

write.table(foiePositive_adjusted,"foiePositiveAdjusted.txt")
