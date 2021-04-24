
### Set working directory

setwd("c:/users/david/dropbox/sad2020/session4/data")

### Required packages

install.packages("nnet") # Only if needed
require(nnet)            # For multinomial logistic regression

install.packages("questionr") # Only if needed
require(questionr)            # For odds-ratio

install.packages("aod") # Only if needed
require(aod)            # For splitbin

install.packages("car") # Only if needed
require(car)            # For Anova

install.packages("bestglm")   # Only if needed
require(bestglm)              # For Model Selection

install.packages("RcmdrMisc")   # Only if needed
require(RcmdrMisc)              # For Stepwise Model Selection

### Import data

# Import data stored in file maturity.txt
dta = read.table("maturity.txt")
# Convert the Maturity variable into a 3-level factor
dta$Maturity = factor(dta$Maturity,levels=c('1','2','3'))
# Overview of data
str(dta)

### Model comparison

# Fit the multinomial logistic model with interaction
maturity.mlogit1 = multinom(Maturity~a*Variety,data=dta)
D1 = deviance(maturity.mlogit1)
D1

# Fit the multinomial logistic model without interaction
maturity.mlogit0 = multinom(Maturity~a+Variety,data=dta)
D0 = deviance(maturity.mlogit0)
D0

# Likelihood-Ratio Test (LRT)
D = D0-D1
D

# Analysis of deviance table for the interaction effect
anova(maturity.mlogit0,maturity.mlogit1,test="Chisq")

# Type-II Analysis of deviance table
Anova(maturity.mlogit1)

### Exhaustive search of the best model for a two-class response

## 2-group maturity

dta12 = subset(dta, Maturity != "3")
# Keeps only rows with Maturity != "3"
dta12 = droplevels(dta12)
# Forgets the unused levels
str(dta12)

# Fit the model with all explanatory variables except 'Variety'
maturity.logit = glm(Maturity~.,family=binomial,data=dta12[,-7])

# Exhaustive search of the best submodel
# Warning Xy = data matrix with y in the last column
maturity.select = bestglm(Xy=dta12[,-7],family=binomial,
                          method="exhaustive",IC="BIC")
maturity.select$Subsets

# Stepwise search of the best submodel
maturity.step = stepwise(maturity.logit,direction="forward/backward",criterion="AIC")
summary(maturity.step)$coefficients

## 3-group maturity

# Fit the model with all explanatory variables except 'Variety'
maturity.mlogit = multinom(Maturity~.,data=dta[,-7])

# Stepwise search of the best submodel
maturity.step = stepwise(maturity.mlogit,direction="forward/backward",criterion="AIC")
summary(maturity.step)


