
## Décider à partir de données  

### Importation des données ELISA
  
elisa <- read.csv2(file="data/elisa_data.csv")
# Structure of object 'elisa'
str(elisa)

# Prints first 6 rows
head(elisa)

# Prints first 8 rows
head(elisa,n=8)

### Premières manipulations de données dans R
  
# Shows the value in row 1, column 2
elisa[1,2] 
  
# Shows all values in column 2
elisa[,2] 
  
# Shows all values in column named "Potency"
elisa[,"Potency"] 
  
# Shows all values in column named "Potency" (same as above)
elisa$Potency 
  
# Shows first 10 values in column named "Potency"
elisa[1:10,"Potency"] 
  
# Shows values for Analyst "1"
elisa[elisa$Analyst==1,"Potency"] 
  
# Data for Analyst "1"
elisa_1 <- elisa[elisa$Analyst=="1",]
elisa_1 <- subset(elisa,Analyst=="1")

### Assign proper data types to variables
  
elisa$Analyst <- factor(elisa$Analyst)
elisa$Run <- factor(elisa$Run)
str(elisa)

### Data summaries
  
# Counts in each group
table(elisa$Analyst,elisa$Run)
  
# Summary statistics for object elisa_1
summary(elisa_1)
  
# Mean, standard deviation (s.d) and relative s.d of potency values
mean(elisa_1$Potency)
sd(elisa_1$Potency)
rel_sd <- sd(elisa_1$Potency)/mean(abs(elisa_1$Potency))*100
  
# Another summary function in package RcmdrMisc
RcmdrMisc::numSummary(data=elisa_1$Potency,
                        statistics=c("mean","sd","quantiles"),
                        groups=elisa_1$Run) 

### Exploration de données, premières synthèses graphiques 
  
# Box plot
boxplot(x=elisa_1$Potency,
          main="Analyst 1",
          ylab="Potency",
          col="lightblue")
grid()
  
# Dot plot
plot(elisa_1$Potency,
       pch=16,
       main="Analyst 1",
       xlab="Value number",ylab="Potency")
grid()
  
# Run-specific boxplots
plot(Potency~Run,
       data=elisa_1,
       main="Analyst 1",
       ylab="Potency")
grid()

## Test d'hypothèse : principes généraux
  
### Exemple 1 : le test de Student

elisa_1 <- subset(elisa,Analyst=="1")
elisa_2 <- subset(elisa,Analyst=="2")

t.test(elisa_2$Potency, mu=1)

power.t.test(n=12,delta=0.1,sd=sd(elisa_2$Potency),type="one.sample")

### Exemple 2 : le test de Student de comparaison de moyennes
  
t.test(Potency~Analyst,data=elisa,var.equal=TRUE)

### Le test de Dixon dans R (package outliers)

dixon_data <- read.csv2("C:/Users/causeur/Dropbox/Axplora/dixon_data.csv")
str(dixon_data)

values <- with(data=dixon_data,Values[Group=="Tests_1_IP"])
# Stores values in group "Tests_1_IP"

library(outliers) # Contains the function dixon.test()

dixon.test(values) # Most extreme is expected to be on both sides
dixon.test(values,two.sided=FALSE) # Most extreme can only happen on one side

## Précision intermédiaire d'une méthode analytique
  
### Analyse de la variance (ANOVA) : principes généraux
  
library(VCA) # Load the R package VCA for variance component analysis
# Create a new variable "Group" as the interaction of Analyst and Run
elisa$Group <- with(elisa,interaction(Analyst,Run))
# Create a new variable "logPotency" as the logarithm of Potency
elisa$logPotency <- log(elisa$Potency)
# Plot the boxplot of logPotency by Group
varPlot(form=logPotency~Analyst/Run, Data=elisa)
grid()

### Test de l'hypothèse nulle : les moyennes par analyste sont les mêmes
  
library(lmerTest)

fit <- lmer(logPotency~Analyst+(1|Analyst:Run),data=elisa)
anova(fit)

### Précision intermédiaire : composantes de la variance
  
fit.IP <- fitVCA(logPotency~Group,Data=elisa,method="anova")
fit.IP

### Précision intermédiaire

sigma2_IP <- fit.IP$aov.tab["total","VC"]
CV <- 100*sqrt(exp(sigma2_IP)-1)
GCV <- 100*(exp(sqrt(sigma2_IP))-1)
