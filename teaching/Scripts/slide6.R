---
title: "Démarche statistique"
subtitle: "Session 6 - Sélection de modèles"
author:
- David Causeur\newline
- Institut Agro Rennes Angers\newline
- IRMAR UMR 6625 CNRS
date: "`r format(Sys.time(), '%d %B, %Y')`"
make149: true
output:
  beamer_presentation: 
    slide_level: 3
  slidy_presentation:
    highlight: tango
  pdf_document:
    highlight: tango
    df_print: kable
  ioslides_presentation:
    toc: yes
    toc_depth: 3
    transition: faster
    widescreen: yes
    number_sections: no
    highlight: tango
    df_print: kable
  html_document:
    toc: yes
    toc_depth: '3'
    df_print: paged
header-includes:
- \usepackage{booktabs}
- \usepackage{wrapfig}
- \usepackage{subcaption}
- \usepackage[font=small,labelfont=bf,labelformat=empty]{caption}
- \usepackage{amsmath}
- \usepackage[T1]{fontenc}
- \usepackage{fancyhdr}
- \usepackage{bm}
- \usepackage{tcolorbox}
- \pagestyle{fancy}
- \fancyhead{}
- \fancyfoot{}
- \fancyfoot[R]{\thepage}
- \setbeamertemplate{footline}{\thepage}
- \def\begincols{\begin{columns}}
- \def\endcols{\begin{columns}}
- \def\begincol{\begin{column}}
- \def\endcol{\begin{column}}
linestretch: 0.8
fontsize: 8pt
vignette: null
---

```{r knitr_init, echo = FALSE, cache = FALSE}
## Global options
options(show.signif.stars = FALSE, width = 50)
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE,
                      prompt = FALSE,
                      tidy = FALSE,
                      comment = NA,
                      message = FALSE,
                      warning = FALSE,
                      fig.align = 'center',
                      fig.width = 7,
                      fig.height = 4.33,
                      size = "tiny",
                      fig.path='pics/', 
                      fig.show='asis')
```

```{r functions, include=FALSE}
# A function for captioning and referencing images
fig <- local({
    i <- 0
    ref <- list()
    list(
        cap=function(refName, text) {
            i <<- i + 1
            ref[[refName]] <<- i
            paste("Figure ", i, ": ", text, sep="")
        },
        ref=function(refName) {
            ref[[refName]]
        })
})
``` 

```{r, echo=FALSE, message=FALSE, warning=FALSE}
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
```

## Modèle linéaire

### Une grande diversité de problématiques

\noindent \textcolor{orange}{Illustration} : 

* La valeur commerciale d'une carcasse de porc est définie à partir de son taux de viande maigre (LMP)
* Mesurer le LMP d'une carcasse de porc est très coûteux (plusieurs heures de dissection)
* On cherche donc à l'approcher au mieux à partir d'informations plus accessibles (épaisseurs de tissus gras et maigres)
* Pour construire une règle de prédiction du LMP, on dispose de mesures sur 354 carcasses

\bigskip

\small

```{r import,echo=TRUE}
dta <- read.table("./data/pig10.txt",stringsAsFactors=TRUE)
str(dta)
```
\normalsize

### Effet groupe

\noindent \textcolor{orange}{Illustration} : les teneurs en viande maigre moyennes par génotype sont-elles les mêmes ?

\begin{columns}
\begin{column}{0.5\textwidth}

\noindent \textcolor{orange}{Analyse de la variance à un facteur}

\bigskip

Pour $i=1,2,3$, $j=1,\ldots,n_{i}$,
\begin{eqnarray*}
Y_{ij} & = & \mu + \alpha_{i} + \varepsilon_{ij} ,
\end{eqnarray*}

où

\begin{itemize}
\item $Y_{ij}$ : teneur en viande maigre de la $j$ème carcasse de génotype $i$
\item $\varepsilon_{ij} \sim {\cal N} ( 0 ; \sigma )$
\item $\alpha_{1}=0$: $\mu$ est la teneur en viande maigre moyenne des carcasses de génotype P$_{0}$ 
\end{itemize}

\end{column}
\begin{column}{0.5\textwidth}  %%<--- here

\hfill
```{r,out.width=150,fig.asp=1,echo=FALSE}
plot(LMP~GENOTYPE,data=dta,
     xlab="Génotype",ylab="LMP",col="orange")
```
\end{column}
\end{columns}

### Analyse de la variance à un facteur dans R

\noindent \textcolor{orange}{Ajustement du modèle}

\medskip

\small

```{r anova1way,echo=TRUE}
mod <- lm(LMP~GENOTYPE,data=dta)
```

\normalsize

\noindent \textcolor{orange}{Test de Fisher} : 2 options (équivalentes)

\medskip

\small

```{r Fopt1,echo=TRUE}
anova(mod)
```
\normalsize

\small

\medskip

```{r Fopt2,echo=TRUE}
mod0 <- lm(LMP~1,data=dta)
anova(mod0,mod)
```
\normalsize

### Description d'un effet groupe

\noindent \textcolor{orange}{Tests post-hoc} : tests de Student de comparaison des génotypes par paires

\medskip

\small

```{r,out.width=180,fig.asp=1,echo=TRUE}
posthoc <- meansComp(mod, ~ GENOTYPE,adjust="bonferroni",graph=TRUE)
```

\normalsize

### Effet d'interaction entre deux variables catégorielles

\noindent \textcolor{orange}{Illustration} : les différences entre les teneurs en viande maigre moyennes par génotype sont-elles les mêmes pour les mâles et les femelles ?

\begin{columns}
\begin{column}{0.5\textwidth}

\noindent \textcolor{orange}{Analyse de la variance à deux facteurs}

\bigskip

Pour $i=1,2,3$, $j=1,2$, $k=1,\ldots,n_{ij}$,
\begin{eqnarray*}
Y_{ijk} & = & \mu + \alpha_{i} + \beta_{j} + (\alpha\beta)_{ij} + \varepsilon_{ijk} ,
\end{eqnarray*}

où

\begin{itemize}
\item $Y_{ijk}$ : teneur en viande maigre de la $k$ème carcasse ayant le génotype $i$
et de sexe $j$
\item $\varepsilon_{ijk} \sim {\cal N} ( 0 ; \sigma )$
\item $\alpha_{1}=0$, $\beta_{1}=0$
\item $(\alpha\beta)_{1j}=0$ pour $j=1,2$, $(\alpha\beta)_{i1}=0$ pour $i=1,2,3$
\end{itemize}

\end{column}
\begin{column}{0.5\textwidth}  %%<--- here

\hfill
```{r,out.width=150,fig.asp=1,echo=FALSE}
interaction.plot(x.factor=dta$GENOTYPE,
                 trace.factor=dta$SEX,bty="l",
                 trace.label="SEXE",type="b",pch=16,
                 response=dta$LMP,lwd=2,lty=1,
                 xlab="Génotype",ylab="LMP",
                 col=c("orange","blue"))
grid()
```

\end{column}
\end{columns}

### Analyse de la variance à deux facteurs dans R

\noindent \textcolor{orange}{Ajustement du modèle avec interaction}

\medskip

\small

```{r,echo=TRUE}
mod <- lm(LMP~SEX+GENOTYPE+GENOTYPE:SEX,data=dta)
anova(mod)
```

\normalsize

\noindent \textcolor{orange}{Ajustement du modèle sans interaction}

\medskip

\small

```{r,echo=TRUE}
mod <- lm(LMP~SEX+GENOTYPE,data=dta)
anova(mod)
```

### Description d'un effet groupe par les moyennes ajustées

\noindent \textcolor{orange}{Tests post-hoc} : tests de Student de comparaison des génotypes par paires

\medskip

\small

```{r,out.width=180,fig.asp=1,echo=TRUE}
posthoc <- meansComp(mod, ~ GENOTYPE,adjust="bonferroni",graph=TRUE)
```

\normalsize


### Effet linéaire

\noindent \textcolor{orange}{Illustration} : la teneur en viande maigre d'une carcasse dépend-elle de son épaisseur de gras ?

\begin{columns}
\begin{column}{0.5\textwidth}

\noindent \textcolor{orange}{Régression linéaire simple}

\bigskip

Pour $i=1,\ldots,n$,
\begin{eqnarray*}
Y_{i} & = & \beta_{0} + \beta_{1} x_{i} + \varepsilon_{i} ,
\end{eqnarray*}

où

\begin{itemize}
\item $Y_{i}$ : teneur en viande maigre de la $i$ème carcasse
\item $x_{i}$ : épaisseur de gras de la $i$ème carcasse
\item $\varepsilon_{i} \sim {\cal N} ( 0 ; \sigma )$
\end{itemize}

\end{column}
\begin{column}{0.5\textwidth}  %%<--- here

\hfill
```{r,out.width=150,fig.asp=1,echo=FALSE}
plot(LMP~LR23Fat,data=dta,pch=16,bty="l",cex.lab=1.25,
     xlab="Epaisseur de gras (mm) - LR23",ylab="LMP",col="orange")
grid()
```
\end{column}
\end{columns}

### Régression linéaire dans R

\noindent \textcolor{orange}{Ajustement du modèle}

\medskip

\small

```{r,echo=TRUE}
mod <- lm(LMP~LR23Fat,data=dta)
```

\normalsize

\noindent \textcolor{orange}{Test de Fisher} : 2 options (équivalentes)

\medskip

\small

```{r,echo=TRUE}
anova(mod)
```
\normalsize

\small

\medskip

```{r,echo=TRUE}
mod0 <- lm(LMP~1,data=dta)
anova(mod0,mod)
```
\normalsize

### Régression linéaire par groupes

\noindent \textcolor{orange}{Illustration} : la relation entre la teneur en viande maigre d'une carcasse et son épaisseur de gras est-elle la même pour les mâles et les femelles ?

\begin{columns}
\begin{column}{0.5\textwidth}

\noindent \textcolor{orange}{Modèle linéaire avec interaction entre l'épaisseur de gras et le sexe}

\bigskip

Pour $i=1,2$, $j=1,\ldots,n_{i}$,
\begin{eqnarray*}
Y_{ij} & = & \beta_{0} + \alpha_{i} + ( \beta_{1} + \gamma_{i} ) x_{ij} + \varepsilon_{ij} ,
\end{eqnarray*}

où

\begin{itemize}
\item $Y_{ij}$ : teneur en viande maigre de la $j$ème carcasse ayant le sexe $i$
\item $x_{ij}$ : épaisseur de gras de la $j$ème carcasse ayant le sexe $i$
\item $\varepsilon_{ij} \sim {\cal N} ( 0 ; \sigma )$
\item $\alpha_{1}=0$, $\gamma_{1}=0$
\end{itemize}

\end{column}
\begin{column}{0.5\textwidth}  %%<--- here

\hfill
```{r,out.width=150,fig.asp=1,echo=FALSE}
dta %>% ggplot() + 
  aes(x = LR23Fat, y = LMP, col = SEX) +
  geom_point() +
  xlab('Epaisseur de gras (mm) - LR23') +
  ylab('LMP') + 
  geom_smooth(method = 'lm')
```

\end{column}
\end{columns}

### Régression linéaire par groupes dans R

\noindent \textcolor{orange}{Ajustement du modèle}

\medskip

\small

```{r,echo=TRUE}
mod <- lm(LMP~LR23Fat+SEX+LR23Fat:SEX,data=dta)
```

\normalsize

\noindent \textcolor{orange}{Test de Fisher} : 

\medskip

\small

```{r,echo=TRUE}
options(width = 300)
anova(mod)
```

\normalsize

## Choix de modèles

### Choix entre deux variables explicatives

\noindent \textcolor{orange}{Illustration} : quelle épaisseur de gras pour expliquer les variations de la teneur en viande maigre, LR23 ou LR34 ? 
\medskip

\noindent \textcolor{orange}{Modèles avec une variable explicative}

\medskip

\small

```{r}
fat_1 <- lm(LMP~LR34Fat,data=dta)
anova(fat_1)
fat_2 <- lm(LMP~LR23Fat,data=dta)
anova(fat_2)
```

\normalsize

\noindent \textcolor{orange}{Conclusion} : l'épaisseur de gras mesurée en LR23 explique mieux les variations de LMP que celle mesurée en LR34 

### Choisir entre deux variables explicatives ou garder les deux ?

\noindent \textcolor{orange}{Modèle avec deux variables explicatives}

\medskip

\small

```{r}
fat_12 <- lm(LMP~LR34Fat+LR23Fat,data=dta)
anova(fat_12)
```

\medskip

```{r}
Anova(fat_12)
```   
   
\normalsize   

\noindent \textcolor{orange}{Deux conclusions antagonistes} 

* Le modèle est mieux ajusté (sommes des carrés des résidus plus faible)
* Une seule des deux épaisseurs de gras est suffisante (LR23) pour expliquer les variations de la teneur en viande maigre 

### Choix de la meilleure variable explicative

\noindent \textcolor{orange}{Illustration} : quelle épaisseur de gras ou de muscle explique le mieux les variations de la teneur en viande maigre ?
   
\bigskip
   
\small
   
```{r,out.width="60%"}
R2 <- cor(dta$LMP,dta[,-c(1:2,10)])^2
ord <- order(R2,decreasing=TRUE)
barplot(R2[ord],col="orange",names.arg = colnames(R2)[ord],las=3,
        ylab=expression(R^2),yaxp=c(0,0.8,8),
        main="Sélection de la meilleure variable explicative")
grid()
```
\normalsize

### Choix du meilleur ensemble de k variables explicatives

   \noindent \textcolor{orange}{Quel modèle ${\cal M}_{k}$ avec $k \leq K$ variables explicatives est suffisant pour expliquer les variations de teneurs en viande maigre ?}
   
   \medskip

   \noindent Choix parmi $2^{K}$ modèles (ici, $K=7$, soit 128 modèles)
   
   \medskip
   
\small

```{r,out.width="60%"}
best <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7)
R2 <- summary(best)$rsq
barplot(R2,col="orange",xlab="Nombre de variables explicatives",
        names.arg = 1:7,ylab=expression(R^2),yaxp=c(0,1,10),
        main="Sélection du meilleur ensemble de k variables explicatives")
grid()
```

\normalsize

### Meilleurs modèles à k variables explicatives

\noindent \textcolor{orange}{Liste des meilleurs modèles à $k$ variables explicatives, $k=1,\ldots,7$}

\medskip

\small

```{r}
options(width=300)
summary(best)$which
```

\normalsize

\medskip

\noindent \textcolor{orange}{Comment choisir parmi ces 7 modèles ? Combien de variables explicatives doit-on garder ?}

### Compromis entre qualité d'ajustement et complexité du modèle

\noindent \textcolor{orange}{Modèle linéaire à $k$ variables explicatives et $p_{k}$ paramètres}

\begin{eqnarray*}
& & \text{Pour } i=1,\ldots,n, \ Y_{i} = \beta_{0} + \beta_{1} x_{i1} + \ldots + \beta_{k} x_{ik} + \varepsilon_{i} , \text{où }
\varepsilon_{i} \sim {\cal N} ( 0 ; \sigma )
\end{eqnarray*}

\medskip

\noindent \textcolor{orange}{Critères d'information}

\medskip

* Akaike Information Criterion (AIC), pour prédire

\begin{eqnarray*}
\text{AIC} & = & n \text{log} \Bigl( \frac{RSS}{n} \Bigr) + 2 p_{k}
\end{eqnarray*}

* Bayesian Information Criterion (BIC), pour expliquer

\begin{eqnarray*}
\text{BIC} & = & n \text{log} \Bigl( \frac{RSS}{n} \Bigr) + \text{log} ( n ) p_{k}
\end{eqnarray*}

\medskip

\noindent \textcolor{orange}{Remarque} : comme $\text{log} ( n ) > 2$ en général, le modèle ayant le BIC le plus faible contient moins de variables que celui ayant l'AIC le plus faible 

### Optimisation des critères d'information dans R

\small

```{r,out.width="75%"}
bic <- summary(best)$bic
aic <- bic+(2-log(nrow(dta)))*(2:8)
plot(1:7,bic,pch=16,bty="l",lwd=2,type="b",col="orange",ylab="IC",
     ylim=range(c(aic,bic)),xlab="Nombre de variables explicatives")
points(1:7,aic,type="b",pch=16,col="blue")
legend("topright",bty="n",col=c("orange","blue"),lwd=2,
       legend=c("BIC","AIC"))
grid()
```

\normalsize

\textcolor{orange}{Remarque} : la 2ème meilleure variable explicative (LR34Fat) n'est pas choisie

### Modèles optimaux

\noindent \textcolor{orange}{Modèle ayant le BIC le plus faible}

\medskip

\small

```{r}
best_bic <- lm(LMP~SplitFat+SplitMuscle+LV23Fat+LR23Fat,data=dta)
Anova(best_bic)
```

\normalsize

\noindent \textcolor{orange}{Modèle ayant le BIC le plus faible}

\medskip

\small

```{r}
best_aic <- lm(LMP~SplitFat+SplitMuscle+LV23Fat+LR23Fat+LR23Muscle,data=dta)
Anova(best_aic)
```

\normalsize

### Méthodes alternatives de recherche du meilleur modèle

\noindent \textcolor{orange}{Méthodes dites pas-à-pas (stepwise) si $K>50$}

   \medskip
   
   \noindent \textbf{Recherche ascendante}
   
   \medskip

   \begin{itemize}
   \item \textbf{Etape 1} : ${\cal M}^{\star}_{1}$, meilleur modèle à une variable explicative

   \medskip

   \item \textbf{Etape $k$} : ${\cal M}^{\star}_{k}$, meilleur modèle parmi ceux complétant ${\cal M}^{\star}_{k-1}$ en ajoutant une variable explicative.

   \medskip

   \item \textbf{Stop} si le BIC de ${\cal M}^{\star}_{k}$ est plus grand que celui de ${\cal M}^{\star}_{k-1}$ .
   
   \end{itemize}
   
   \medskip
   
   \noindent \textbf{Recherche descendante}
   
   \medskip

   \begin{itemize}
   \item \textbf{Etape 1} : ${\cal M}^{\star}_{K}$, modèle contenant toutes les variables explicatives

   \medskip

   \item \textbf{Etape $k$} : ${\cal M}^{\star}_{K-k+1}$, meilleur modèle parmi ceux obtenus à partir de ${\cal M}^{\star}_{K-k+2}$ en enlevant une variable explicative.

   \medskip

   \item \textbf{Stop} si le BIC de ${\cal M}^{\star}_{K-k+1}$ est plus grand que celui de ${\cal M}^{\star}_{K-k+2}$ .
   
   \end{itemize}

### Méthodes pas-à-pas dans R

\noindent \textcolor{orange}{Recherches exhaustive, ascendante et descendante}

\medskip

\small

```{r}
best_exh <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="exhaustive")
best_fwd <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="forward")
best_bwd <- regsubsets(LMP~.,data=dta[,-(1:2)],nvmax=7,
                   method="backward")
```

\normalsize

\medskip

\noindent \textcolor{orange}{Remarques} : 

* Dans le cas présent, les trois options de recherche conduisent au même résultat
* Toutefois, en général, la recherche pas-à-pas ne garantit pas de trouver le meilleur modèle

## Ce qu'il faut retenir

### Messages principaux

\noindent \textcolor{orange}{Lorsque la problématique cible précisément l'effet d'une variable explicative}

* Recenser les confusions potentielles avec l'effet de cette variable explicative et les effets d'interaction avec d'autres variables explicatives
* Procéder par des tests de Fisher de comparaison de modèles 

\medskip

\noindent \textcolor{orange}{Lorsque la problématique ne cible pas en particulier une variable explicative}

* Comparer les modèles construits à partir de tous les sous-ensembles possibles de variables explicatives
* Choisir le modèle optimisant un critère d'information (AIC ou BIC) pour réaliser le meilleur compromis entre qualité d'ajustement et complexité du modèle

\medskip

\noindent \textcolor{orange}{Pour aller plus loin} : voir [le principe philosophique du rasoir d'Ockham](https://en.wikipedia.org/wiki/Occam%27s_razor)
