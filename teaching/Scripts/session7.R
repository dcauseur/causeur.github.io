
Install_and_Load <- function(packages) {
  k <- packages[!(packages %in% installed.packages()[, "Package"])];
  if(length(k)) {
    install.packages(k, repos = "https://cran.rstudio.com/");
  }
  for(package_name in packages) {
    library(package_name, character.only = TRUE, quietly = TRUE);
  }
}

Install_and_Load(c("kableExtra","ggpubr","car",
                   "RcmdrMisc","FactoMineR","ggplot2",
                   "wesanderson","tidyverse","dplyr","GGally",
                   "ggcorrplot"))

### Objectif : visualiser de manère simplifiée des profils décrits pas plusieurs variables (quantitatives)

dta <- read.table("./data/pig10.txt",stringsAsFactors=TRUE)
dta %>% select(c(3,5,6,8,4,7,9)) %>% 
  cor() %>% round(1) %>% ggcorrplot(type = "lower",lab = TRUE)

### Illustration pour deux variables fortement corrélées  

pca_1 <- PCA(dta[,c(6,8)],graph=FALSE)
sc_x <- scale(dta[,3:9]) 
plot(LR23Fat~LR34Fat,data=sc_x,bty="n",pch=16,pty="s",
     col="orange",xlab="Epaisseur de gras (mm) - LR34",
     ylab="Epaisseur de gras (mm) - LR23",ylim=c(-3.75,3.75),xlim=c(-3.75,3.75))
abline(h=0,v=0,col="darkgray",lwd=2)
abline(a=0,b=1,col="blue",lwd=2)
text(-3,-2.6,expression(z[1]),cex=1,col="blue",adj=0)
text(-3,2.6,expression(var(z[1])==1.95),lwd=2,cex=1,adj=0)
grid()

### Illustration pour deux variables modérément corrélées  

pca_1 <- PCA(dta[,c(4,9)],graph=FALSE)
sc_x <- scale(dta[,3:9]) 
plot(LR34Muscle~SplitMuscle,data=sc_x,bty="n",pch=16,pty="s",
     col="orange",xlab="Epaisseur de muscle (mm) - Split",
     ylab="Epaisseur de muscle (mm) - LR34",ylim=c(-3.75,3.75),xlim=c(-3.75,3.75))
abline(h=0,v=0,col="darkgray",lwd=2)
abline(a=0,b=1,col="blue",lwd=2)
text(-3,-2.6,expression(z[1]),cex=1,col="blue",adj=0)
text(-3,2.6,expression(var(z[1])==1.46),lwd=2,cex=1,adj=0)
grid()

### Inertie d'une composante principale  

pca <- PCA(dta[,c("LR23Fat","LR34Fat")],graph=FALSE)
round(pca$eig,3)
### Analyse en Composantes Principales (ACP) de deux variables

plot(pca,choix="ind",pch=16,col.ind="orange",label="none")

### Interprétation des composantes principales par les scores extrêmes

z1 <- pca$ind$coord[,1]
select_high <- which(z1>4)
select_low <- which(z1< -2.5)
cbind(dta[,c(1,2,6,8,10)],PC1=z1)[c(select_high,select_low),]

### Interprétation des composantes principales par le cercle des corrélations

plot(pca,choix="var")

### Analyse en Composantes Principales (ACP) de p variables

pca <- PCA(dta[,3:9],graph=FALSE)
barplot(pca$eig[,3],col="orange",main="Inertie expliquée cumulée (%)")
grid()

### Plan factoriel 

plot(pca,choix="ind",pch=16,col.ind="orange",label="none",axes=c(1,2))

### Qualité de représentation d'un individu dans le plan factoriel 

outlier <- which.max(pca$ind$coord[,1])
dta[outlier,-(1:2)]
pca$ind$cos2[outlier,]

### Interprétation des composantes principales 

plot(pca,choix="var")

### Qualité de représentation d'une variable dans le plan factoriel 

round(pca$var$coord["SplitMuscle",],digits=3)

pca$var$cos2["SplitMuscle",]

### Interprétation à l'aide de variables catégorielles supplémentaires 

pca <- PCA(dta,quali.sup=1:2,quanti.sup=10,graph=FALSE)
plot(pca,choix="ind",label="none",col.hab=c("orange","blue"),habillage=2)

round(pca$quali.sup$eta2,3)

### Interprétation à l'aide de variables quantitatives supplémentaires 

plot(pca,choix="var")

round(pca$quanti.sup$cor,3)

