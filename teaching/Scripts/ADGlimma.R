
# Installation of limma (only if needed)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma",update=FALSE)

require(limma)

### Set working directory (where data are)

setwd("C:/Users/David/Dropbox/adg2020/Slides/Session2/Data")
dir()   # Displays the files in directory

### Experimental design

targets = readTargets("targetsF.txt")
head(targets)    # Display the first 6 rows

dietgenotype = factor(targets$Factor)
table(dietgenotype)

genotype = factor(substring(dietgenotype,first=1,last=1))
diet = factor(substring(dietgenotype,first=3,last=4))

table(diet,genotype)

### Import external data

# The Excel file on the platform has been converted into csv (seperator=;) 
external_data = read.csv2("44variables_poulet9S_V1.csv")
str(external_data)

# The first column contains the number of the microarray
# It should be consistent with the number in the ArrayName column of targets

# Extracts the microarray number from ArrayName in targets
ArrayNumber = substring(targets$ArrayName,first=3,last=4)
ArrayNumber = as.numeric(ArrayNumber)

# Checks the consistency with the numbers in the 1st column of external_data
all(ArrayNumber==external_data[,1])

### Import data

# Modify the name of the files to be imported (remove the path)
arrayfile = targets$ArrayFile
targets$ArrayFile = substring(arrayfile,first=13,last=nchar(arrayfile))   

# Quality control 

MyFiltering <- function(x)   {
  okType = x$ControlType==0          # Probe=gene
  okExpress = x$gIsWellAboveBG==1    # Gene expression >> background
  # No flag + no spatial aggregates
  okSpotQuality = (x$IsManualFlag ==0)&(x$gIsFeatNonUnifOL == 0) 
  as.numeric(okType & okExpress & okSpotQuality)
}

# Import data

G = read.maimages(files = targets$ArrayFile, names=targets$ArrayName, 
                  source="agilent",green.only=TRUE,wt.fun=MyFiltering, 
                  columns=list(E = "gMedianSignal", Eb = "gBGMedianSignal"))

typeof(G)   # G is a list 		
names(G)		# G has 6 components 

str(G$genes) # G$genes gives information on probes    
dim(G$E)     # G$E is the probe x microarray matrix of gene expressions 
head(G$E[,1:10]) # Display the first 6 row and 10 columns

# Restricting data to gene expressions

G = G[G$genes$ControlType==0,]
dim(G) 

### Data quality

dim(G$weights)    
head(G$weights[,1:10])  # G$weights is a probe x microarray matrix of boolean values
                        # 0 = bad spot, 1 = good spot    

arrayquality = colMeans(G$weights) # Proportions of good spots on each microarray 
summary(arrayquality)

barplot(arrayquality,ylim=c(0,1), las=3,main="% de spots de bonne qualité") 
abline(h=0.7, col="red", lwd=2)

out = which(arrayquality>0.9)     # Suspicious microarrays
G = G[,-out]                      # Two microarays are removed
dietgenotype = dietgenotype[-out]
diet = diet[-out]
genotype = genotype[-out]
targets = targets[-out,]
external_data = external_data[-out,]

# Proportions of good spots for each probe
genesqualityGBL = rowMeans(G$weights[,dietgenotype=="G_BL"])
genesqualityMBL = rowMeans(G$weights[,dietgenotype=="M_BL"])
genesqualityGHL = rowMeans(G$weights[,dietgenotype=="G_HL"])
genesqualityMHL = rowMeans(G$weights[,dietgenotype=="M_HL"])

# Good probe if at least 75% of good spots in at least one experimental plot 
okgenes = (genesqualityGBL>=0.75)|
          (genesqualityMBL>=0.75)|
          (genesqualityGHL>=0.75)|
          (genesqualityMHL>=0.75)

sum(okgenes)

# Restriction to good probes

G = G[okgenes,]
dim(G)           # 46741 probes x 46 microarrays

### Data normalization

# log2 transformation

G$E = log2(G$E)     # log-transformation of expression data
G$Eb = log2(G$Eb)   # log-transformation of background values

# Examination of background

boxplot(G$Eb,main="Background",names=targets$ArrayName,las=3)

G = backgroundCorrect(G, method="none") # Another possibility is method="subtract"
names(G) # The background component has been removed. The matrix E is unchanged.

# Examination of gene expressions

boxplot(G$E,main="Signal",names=targets$ArrayName,las=3)

# Median normalization

medians = apply(G$E,2,median)  # Medians for each microarray 
G$E = sweep(G$E,2,medians)     # Subtract medians for each microarray (columns)                               

boxplot(G$E,main="Signal",names=targets$ArrayName,las=3)

# Handling of replicates

head(G$genes$ProbeName)                  # First 6 probe names
names_counts = table(G$genes$ProbeName)  # Numbers of probes for each probe name
table(names_counts)                      # 681 probe names are replicated 

G$E = avereps(G$E, ID=G$genes$ProbeName) # Replace replicated probes by average 

dim(G$E)                          # Dimensions of the normalized expression data                                 

# Export normalized gene expression data

write.table(G$E,"foie.txt")

# Export design information

covariates = data.frame(Diet=diet,Genotype=genotype)
rownames(covariates) = targets$ArrayName
str(covariates)

write.table(covariates,"covariates.txt")

# Export external data

write.table(external_data,"external.txt")
