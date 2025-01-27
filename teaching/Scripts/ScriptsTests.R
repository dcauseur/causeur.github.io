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
                   "wesanderson","tidyverse","dplyr","GGally",
                   "ggcorrplot","broom","boot"))

### Importer des données

fruit <- read.table(file="data/fruit.txt",stringsAsFactors=TRUE)
str(fruit)

### Les poids moyens des fruits par variété sont-il différents ?

tab <- fruit %>% group_by(Variete) %>%
  summarise(Mean=mean(Poids),Sd=sd(Poids),
              Min=min(Poids),Max=max(Poids))
tab

### Représentation des variétés dans l'échantillon
  
fruit %>% ggplot() + 
  geom_bar(aes(x = Variete),fill="orange",width=0.5) +
  ggtitle("Effectifs par variété") + xlab("Variété") + ylab("Effectifs")

### Répartition des poids des fruits dans l'échantillon

fruit %>% ggplot() + aes(x = Poids) + 
  geom_histogram(fill="orange",bins=15,col="white") +
  ggtitle("Répartition des poids") + xlab("Poids") +
  ylab("Effectifs")

### Répartition des poids par variété
fruit %>% ggplot() + 
  geom_boxplot(aes(x=Variete,y=Poids),fill="orange") +
  ggtitle("Répartition des poids par variété") +
  xlab("Variété") + ylab("Poids")

### Modèle d'analyse de la variance à un facteur

mod1 <- lm(Poids~Variete,data=fruit)
coef(mod1)

mod0 <- lm(Poids~1,data=fruit)
coef(mod0)

### Test de Fisher
  
anova(mod0,mod1)

### Comparaison de deux moyennes

fruit2 <- fruit %>% filter(Variete%in%c("37","go")) %>% droplevels()
t.test(Poids~Variete,data=fruit2,var.equal=TRUE)

### Puissance du dispositif experimental

mod = lm(Poids~Variete,data=fruit2)
sigma = summary(mod)$sigma
power.t.test(delta=2,n=100,sd=sigma,sig.level=0.05)$power

### Comparaisons par paires

comp <- emmeans(mod1, ~ Variete)
pairs(comp,adjust="bonf")

### Comparaison graphique des variétés de fruits selon leur poids

res <- meansComp(mod1, ~ Variete,graph=TRUE)

### Diamètre d'un fruit en fonction de son indice de clarté

fruit %>% ggplot() + aes(x=L,y=Diam) + 
  geom_point() + geom_smooth(method='lm') +
  ggtitle("Relation entre diamètre et clarté") + 
  xlab("L") + ylab("Diamètre")

### Corrélation

cor(fruit$L,fruit$Diam)

### Modèle de régression linéaire

mod1 <- lm(Diam~L,data=fruit)
coef(mod1)

summary(mod1)$r.squared

### Evaluation graphique de la qualité de l'ajustement 
  
fruit$Diam_fit <- fitted(mod1)
fruit %>% ggplot(aes(x=Diam,y=Diam_fit)) + geom_point() + geom_smooth(method='lm') +
  ggtitle("Diamètres ajustés et observés") + geom_abline(intercept=0,slope=1) +
  xlab("Diamètres observés") + ylab("Diamètres ajustés") + ylim(34,54)

### Test de Fisher
  
mod0 <- lm(Diam~1,data=fruit)
anova(mod0,mod1)

## Effet linéaire par groupes
  
### Visualisation d'un effet linéaire par groupes
  
p <- fruit %>% ggplot() + geom_point(aes(x=L,y=Diam,col=Variete)) +
    ggtitle("Relation entre diamètre et clarté") + xlab("Clarté") + ylab("Diamètre")
p

### Tendance linéaire par groupes
  
p + geom_smooth(method="lm",aes(x=L,y=Diam,col=Variete))

### Modèle pour un effet linéaire par groupes
  
mod2 <- lm(Diam~L*Variete,data=fruit)
summary(mod2)$coefficients

### Equation unique ou équations par groupe ?

anova(mod1,mod2)

## Choix du meilleur modèle

### Prédiction du taux de sucre d'un fruit
    
mod <- lm(Saccharose~Poids+Diam+L+a+b,data=fruit)
coef(mod)
summary(mod)$r.squared

mod0 = lm(Saccharose~1,data=fruit)
anova(mod0,mod)

### Confusion d'effets
      
car::Anova(mod)

### Choix du meilleur modele

bestmod <- regsubsets(Saccharose~Poids+Diam+L+a+b,data=fruit)
tidy(bestmod) %>% dplyr::select(1:7) 

### Compromis entre qualité d'ajustement et complexité du modèle  
        
tidy(bestmod) %>% dplyr::select(c(1:6,9))

### Sélection pas à pas 
        
mod <- lm(Saccharose~Poids+Diam+L+a+b,data=fruit)
select <- stepwise(mod,direction="forward/backward",
                           criterion="BIC",trace=0)
coef(select)

