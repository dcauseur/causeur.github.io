
### Working directory (where data files are stored)

setwd("C:/Users/David/dropbox/adg2020")

### Import data

chicken = read.table("chicken.txt")

### Data summary

dim(chicken)            # 27 microarrays - 314 genes
str(chicken[,1:10])     # Overview of first 10 gene expressions

### Missing data

na.chicken = is.na(chicken)         # 0/1 matrix locating missing data  
dim(na.chicken)                     # Numbers of rows and columns 
head(na.chicken[,1:10])             # First 6 rows, first 10 variables
na.counts = colSums(isna.chicken)   # Count missing data per gene 

table(na.counts)                    # Summarize the distribution of NA counts

out = na.counts>=5                  # Genes with more than 5 NAs will be excluded
sum(out)                            # How many excluded genes?
chicken = chicken[,-which(out)]     # Exlusion of genes with too many NAs
dim(chicken)                        # Numbers of rows and columns

### Design

design_labels = rownames(chicken)      # Two way design explicit in the rownames
design_labels

# 1st character gives the diet
diet_labels = substring(design_labels,first=1,last=1) 
diet_labels

# 2nd character gives the phenotype
pheno_labels = substring(design_labels,first=2,last=2) 
pheno_labels

# Displays the two-way (almost) balanced design
table(diet_labels,pheno_labels) 

# Create a dataset with design information 
diet = factor(diet_labels)       # convert diet_labels as a factor 
summary(diet_labels)
summary(diet)

phenotype = factor(pheno_labels) # convert phenotype_labels as a factor
covariates = data.frame(dief=diet,phenotype=phenotype)
str(covariates) # covariates has two columns, diet and phenotype

### Test for the diet effect

## t-test for one gene (the first one, arbitrarily)

t.test(y~x,data=data.frame(x=diet,y=chicken[,1]),var.equal=TRUE)

## Equivalent F-test 

mod = lm(y~x,data=data.frame(x=diet,y=chicken[,1]))
anova(mod)       # Analysis of Variance Table
anova(mod)[1,5]  # A focus on the p-value

## Simultaneous tests for all genes

# Create a function giving the p-value of a one-way anova
MyAnova = function(group,gene) {
   mod = lm(gene~group)
   anovatab = anova(mod)
   pval = anovatab[1,5]
   return(pval)
}

# Illustration of how it works
MyAnova(diet,chicken[,1])
MyAnova(diet,chicken[,2])

# Calculation of p-values for the diet effect on all gene expressions
pval = apply(X=chicken,MARGIN=2,FUN=MyAnova,group=diet)

# Positive genes (significant diet effect on expression)
positive = pval<=0.05
sort(pval[positive])   # p-values for positive genes in ascending order

## Can we be sure these genes are differentially expressed?

# Simulation of a null dataset (no diet effect)
dta0 = matrix(rnorm(27*303),nrow=27,ncol=303)
dim(dta0)          # Numbers of rows and columns
head(dta0[,1:10])  # Display the first 6 rows and first 10 columns 

# Calculation of p-values for the diet effect on all variables
pval0 = apply(X=dta0,MARGIN=2,FUN=MyAnova,group=diet)

# Positive variables (significant diet effect)
positive0 = pval0<=0.05
sort(pval0[positive0])
mean(positive0)        # proportion of false positives

# Histogram of null p-values
hist(pval0,xlab="Null p-values",proba=TRUE,main="Distribution of null p-values",
     col="orange",cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(h=1,lwd=2,col="darkgray")

### Signal-to-noise ratio

## Back to the F-test

mod = lm(y~x,data=data.frame(x=diet,y=chicken[,1]))
anova(mod)       # Analysis of Variance Table

## F-test is a signal-to-noise ratio

anova(mod)[1,3]  # Signal (between diet variation of gene expression)
anova(mod)[2,3]  # Noise (within diet variation of gene expression)

## A closer look at signal and noise 

# Function to extract signal
Signal = function(group,gene) {
   mod = lm(gene~group)
   anovatab = anova(mod)
   return(anovatab[1,3])
}

# Function to extract noise
Noise = function(group,gene) {
   mod = lm(gene~group)
   anovatab = anova(mod)
   return(anovatab[2,3])
}

# Calculation of signal and noise values for all gene expressions
signal = apply(X=chicken,MARGIN=2,FUN=Signal,group=diet)
noise = apply(X=chicken,MARGIN=2,FUN=Noise,group=diet)

# Relationship between noise and signal (log2 scale)
plot(log2(signal),log2(noise),pch=16,xlab="Signal (log2 scale)",ylab="Noise (log2 scale)",
   main="Noise versus signal")

### First normalization operation: log2 transformation 

chicken2 = log2(chicken)

## Impact of log2-transformation of gene expression on signal/noise relationship

signal2 = apply(X=chicken2,MARGIN=2,FUN=Signal,group=diet)
noise2 = apply(X=chicken2,MARGIN=2,FUN=Noise,group=diet)

plot(log2(signal2),log2(noise2),pch=16,xlab="Signal (log2 scale)",ylab="Noise (log2 scale)",
   main="Noise versus signal after log2-transformation")

## Simultaneous tests of the diet effect after log2-transformation

pval2 = apply(X=chicken2,MARGIN=2,FUN=MyAnova,group=diet)

positive2 = pval2<=0.05
sort(pval2[positive2])

### Power of the design to test for a diet effect

## What is the probability of detecting a difference of 1 (log2-scale) between 
## gene expressions in the two diets?

# Illustration for gene 1

sqrt(noise2[1]) # Gives the within-diet standard deviation of gene expression
power.t.test(delta=1,sig=0.05,sd=sqrt(noise2[1]),n=14)

# Now for all genes

vpower = power.t.test(delta=1,sig=0.05,sd=sqrt(noise2),n=14)
summary(vpower$power)

### Noise reduction

# Distribution of gene expressions in each microarrays
boxplot(t(chicken2),col="darkgray",bty="l",xlab="",ylab="Gene expression",
main="Between microarrays gene expression variations",cex.axis=1.25,cex.lab=1.25,
cex.main=1.25,las=3)

### Median normalization

# Calculate median gene expression g=for each microrray
medians = apply(X=chicken2,MARGIN=1,FUN=median,na.rm=TRUE)
medians

# Subtract medians of each row (row=microarray) of data table
chicken3 = sweep(x=chicken2,MARGIN=1,FUN="-",STATS=medians)

# Impact of median normalization on the distribution of gene expression in microarrays

summary(as.numeric(chicken3[1,])) # 1st microarray
summary(as.numeric(chicken3[2,])) # 1st microarray

boxplot(t(chicken3),col="darkgray",bty="l",xlab="",ylab="Gene expression",
main="Between microarrays gene expression variations",cex.axis=1.25,cex.lab=1.25,
cex.main=1.25,las=3)
mtext("After median normalization")

# Impact of median normalization on within-diet standard deviations

noise3 = apply(X=chicken3,MARGIN=2,FUN=Noise,group=diet)
summary(sqrt(noise2)) # Standard deviations before median normalization
summary(sqrt(noise3)) # Standard deviations after median normalization

# Impact of median normalization on between-diet variations

signal3 = apply(X=chicken3,MARGIN=2,FUN=Signal,group=diet)
summary(signal2) # Between-diet variation before median normalization
summary(signal3) # Between-diet variation after median normalization

# Impact of median normalization on significance

pval3 = apply(X=chicken3,MARGIN=2,FUN=MyAnova,group=diet)

positive3 = pval3<=0.05
sort(pval3[positive3])

sum(positive3) # Number of positive genes

### An alternative normalization procedure: quantile normalization

# Installation of limma (only if needed)

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("limma",update=FALSE)

require(limma)

# Quantile normalization

chicken4 = normalizeBetweenArrays(t(chicken2),method="quantile")
dim(chicken4)   # In limma, gene expression datasets are in gene x microarray format

# Transpose the data table
chicken4 = t(chicken4)

# Impact of quantile normalization on distributions of gene expression in microarrays
summary(chicken4[1,]) # 1st microarray
summary(chicken4[2,]) # 1st microarray

boxplot(t(chicken4),col="darkgray",bty="l",xlab="",ylab="Gene expression",
main="Between microarrays gene expression variations",cex.axis=1.25,cex.lab=1.25,
cex.main=1.25,las=3)
mtext("After quantile normalization")

# Impact of quantile normalization on significance

pval4 = apply(X=chicken4,MARGIN=2,FUN=MyAnova,group=diet)

positive4 = pval4<=0.05
sort(pval[positive4])

sum(positive4)

# Impact of normalization on power

noise4 = apply(X=chicken3,MARGIN=2,FUN=Noise,group=diet)
vpower = power.t.test(delta=1,sig=0.05,sd=sqrt(noise4),n=14)

summary(vpower$power)

### Comparison of the two normalization procedures
 
# Install VennDiagram (only if needed)

install.packages("VennDiagram") 
library(VennDiagram)

# Display Venn diagram comparing the normalization procedures

compare.norm = venn.diagram(list(Median=which(positive3),
            Quantiles=which(positive4)),filename=NULL,
            col="orange",fill="orange",
            cat.cex=1.5,cex=1.5,main.cex=1.5,
            main.fontfamily = "sans",fontfamily="sans",cat.fontfamily="sans",
            category.names=c("Median \n normalization","Quantile \n normalization"),
            main="Number of positive genes \n using median and quantile normalization")
grid.newpage()
grid.draw(compare.norm)


