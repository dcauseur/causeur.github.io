
# Installation et chargement de packages

Install_and_Load <- function(packages) {
  k <- packages[!(packages %in% installed.packages()[, "Package"])];
  if(length(k)) {
    install.packages(k, repos = "https://cran.rstudio.com/");
  }
  for(package_name in packages) {
    library(package_name, character.only = TRUE, quietly = TRUE);
  }
}

Install_and_Load(c("kableExtra","emmeans","ggpubr","car",
                   "leaps","RcmdrMisc","FactoMineR","ggplot2",
                   "wesanderson","tidyverse","dplyr"))

## Modèle linéaire

### Importation des doonées

dta <- read.table("./data/pig10.txt",stringsAsFactors=TRUE)
str(dta)

### Effet groupe 

plot(LMP~GENOTYPE,data=dta,
     xlab="Génotype",ylab="LMP",col="orange")

### Analyse de la variance à un facteur

mod <- lm(LMP~GENOTYPE,data=dta)
anova(mod)

mod0 <- lm(LMP~1,data=dta)
anova(mod0,mod)

### Description d'un effet groupe

posthoc <- meansComp(mod, ~ GENOTYPE,adjust="bonferroni",graph=TRUE)

### Analyse de la variance à deux facteurs dans R

mod <- lm(LMP~SEX+GENOTYPE+GENOTYPE:SEX,data=dta)
anova(mod)

mod <- lm(LMP~SEX+GENOTYPE,data=dta)
anova(mod)

### Description d'un effet groupe par les moyennes ajustées

posthoc <- meansComp(mod, ~ GENOTYPE,adjust="bonferroni",graph=TRUE)

### Effet linéaire

plot(LMP~LR23Fat,data=dta,pch=16,bty="l",cex.lab=1.25,
     xlab="Epaisseur de gras (mm) - LR23",ylab="LMP",col="orange")
grid()

### Régression linéaire dans R

mod <- lm(LMP~LR23Fat,data=dta)
anova(mod)

mod0 <- lm(LMP~1,data=dta)
anova(mod0,mod)

### Régression linéaire par groupes

mod <- lm(LMP~LR23Fat+SEX+LR23Fat:SEX,data=dta)
anova(mod)

## Choix de modèles

### Choix entre deux variables explicatives

fat_1 <- lm(LMP~LR34Fat,data=dta)
anova(fat_1)
fat_2 <- lm(LMP~LR23Fat,data=dta)
anova(fat_2)

### Choisir entre deux variables explicatives ou garder les deux ?

fat_12 <- lm(LMP~LR34Fat+LR23Fat,data=dta)
anova(fat_12)

Anova(fat_12)

### Choix de la meilleure variable explicative

R2 <- cor(dta$LMP,dta[,-c(1:2,10)])^2
ord <- order(R2,decreasing=TRUE)
barplot(R2[ord],col="orange",names.arg = colnames(R2)[ord],las=3,
        ylab=expression(R^2),yaxp=c(0,0.8,8),
        main="Sélection de la meilleure variable explicative")
grid()

### Choix du meilleur ensemble de k variables explicatives

best <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7)
R2 <- summary(best)$rsq
barplot(R2,col="orange",xlab="Nombre de variables explicatives",
        names.arg = 1:7,ylab=expression(R^2),yaxp=c(0,1,10),
        main="Sélection du meilleur ensemble de k variables explicatives")
grid()

### Meilleurs modèles à k variables explicatives

summary(best)$which

### Compromis entre qualité d'ajustement et complexité du modèle

bic <- summary(best)$bic
aic <- bic+(2-log(nrow(dta)))*(2:8)
plot(1:7,bic,pch=16,bty="l",lwd=2,type="b",col="orange",ylab="IC",
     ylim=range(c(aic,bic)),xlab="Nombre de variables explicatives")
points(1:7,aic,type="b",pch=16,col="blue")
legend("topright",bty="n",col=c("orange","blue"),lwd=2,
       legend=c("BIC","AIC"))
grid()

### Modèles optimaux

best_bic <- lm(LMP~SplitFat+SplitMuscle+LV23Fat+LR23Fat,data=dta)
Anova(best_bic)

best_aic <- lm(LMP~SplitFat+SplitMuscle+LV23Fat+LR23Fat+LR23Muscle,data=dta)
Anova(best_aic)

### Méthodes alternatives de recherche du meilleur modèle

### Méthodes pas-à-pas dans R

best_exh <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="exhaustive")
best_fwd <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="forward")
best_bwd <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="backward")

