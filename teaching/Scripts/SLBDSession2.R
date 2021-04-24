

# Set working directory

setwd("C:/Users/David/Dropbox/ADB_2021/Cours/Session2")

# Required packages

install.packages("nnet")  # Installation is only needed if the package is missing
install.packages("leaps") 
install.packages("RcmdrMisc")
install.packages("pls")
install.packages("groupdata2")
install.packages("boot")

require(nnet)            # Multinomial logistic regression
require(leaps)           # For regsubsets
require(RcmdrMisc)       # For Stepwise
require(pls)             # For segments
require(groupdata2)      # For fold
require(boot)            # For cv.glm

# Regression with a real-valued response

## Import 'invasive probe' dataset

pig = read.table("./data/invasive.txt")
str(pig)  # Overview of data
dim(pig)  # Number of rows and columns in pig

summary(pig)   # Provides a columnwise summary of the data table

## Exhaustive search of the best model

### Rank the models according to their RSS (using package leaps)
select = summary(regsubsets(LMP~.,data=pig,nvmax=11))

### Boolean selection matrix for best submodels 
select$which

### Best sub-model with one variable
colnames(select$which)[select$which[1,]]

### Best sub-model with two variables
colnames(select$which)[select$which[2,]]

### RSS plot for exhaustive feature selection
plot(1:11,select$rss,type="b",pch=16,xlab="Number of variables in the submodel",
     ylab="Residual sum-of-squares",main="Quality of fit of best submodels",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2,col="darkgray",bty="l")
grid()

### Equivalent R2 plot
plot(1:11,select$rsq,type="b",pch=16,xlab="Number of variables in the submodel",
     ylab=expression(R^2),main="Quality of fit of best submodels",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2,col="darkgray",bty="l")
grid()

## Prediction performance of the best submodels  

### Best submodel with one explanatory variable
mod = glm(LMP~.,data=pig[,select$which[1,]])

### 10-fold CV PRESS 
cvmod = cv.glm(pig[,select$which[1,]],mod,K=10)
cvmod$delta           # MSEP
nrow(pig)*cvmod$delta # PRESS
select$rss[1]         # RSS

press = rep(0,11)     # vector of PRESS for best sub-models
for (j in 1:11) {
   mod = glm(LMP~.,data=pig[,select$which[j,]])
   cvmod = cv.glm(pig[,select$which[j,]],mod,K=10)
   press[j] = nrow(pig)*cvmod$delta[2] 
}

# PRESS plot for exhaustive feature selection
plot(1:11,select$rss,type="b",pch=16,xlab="Number of variables in the submodel",
     ylab="Residual sum-of-squares",main="Quality of fit of best submodels",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2,col="darkgray",bty="l")
lines(1:11,press,type="b",col="blue",pch=15,lwd=2)
legend("topright",lwd=2,pch=c(16,15),legend=c("Internal validation","Cross validation"),
       bty="n",cex=1.25,col=c("darkgray","blue"))
grid()

### BIC - AIC

bic = select$bic                            # BIC
aic = bic - (log(nrow(pig))-2)*(2:12)       # AIC

# BIC-AIC plot for exhaustive feature selection
plot(1:11,bic,pch=16,bty="l",type="b",xlab="Number of explanatory variables",
     ylab="Information criterion",ylim=range(c(aic,bic)),col="darkgray",
     main="Exhaustive model selection",cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2)
lines(1:11,aic,type="b",pch=17,lwd=2,col="coral1")
legend("topleft",lwd=2,lty=1,pch=c(16,17),col=c("darkgray","coral1"),bty="n",cex=1.25,legend=c("BIC","AIC"))
grid()

selected = select$which[which.min(bic),]   # Indices of selected variables (+ LMP)
bestmod = glm(LMP~.,data=pig[,selected])   # Fits the best submodel
coef(bestmod)

## Prediction performance of the best submodel (minimizing AIC)

n = nrow(pig)                   # Sample size
segments = cvsegments(k=10,N=n) # Defines a list of 10 random segments
segments

cvpredictions = rep(0,n)   # Initialize a n-vector of predicted LMP 
for (k in 1:10) {
   train = pig[-segments[[k]],]   # Training dataset
   test = pig[segments[[k]],]     # Test dataset
   select = summary(regsubsets(LMP~.,data=train,nvmax=11))
   bic = select$bic                            # BIC
   selected = select$which[which.min(bic),]   # Indices of selected variables (+ LMP)
   bestmod = glm(LMP~.,data=train[,selected])   # Fits the best submodel
   cvpredictions[segments[[k]]] = predict(bestmod,newdata=test)
}

PRESS = sum((pig$LMP-cvpredictions)^2)

# PRESS plot for exhaustive feature selection
plot(1:11,press,type="b",pch=16,xlab="Number of variables in the submodel",
     ylab="Residual sum-of-squares",main="Quality of fit of best submodels",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2,col="darkgray",bty="l")
abline(h=PRESS,col="blue",pch=15,lwd=2)
legend("topright",lwd=2,pch=c(16,16),legend=c("Best submodels","Best submodel"),
       bty="n",cex=1.25,col=c("darkgray","blue"))
grid()

# Regression with a K-class response

## Import coffee data

coffee = read.table("./data/coffee.txt")
dim(coffee)       # Number of rows and columns in data
str(coffee)       # Overview of data
coffee$Localisation = factor(coffee$Localisation) # Convert 'Localisation' into a factor
summary(coffee)   # Provides a columnwise summary of the data table (8 first columns)

## ML fit the most complete model for 'Localisation'

mod = multinom(Localisation~.,data=coffee,maxit=200,trace=FALSE) 
   # ML fit of the logistic model

## Stepwise search of the best model using package RcmdrMisc

select = stepwise(mod,direction="forward/backward",criterion="AIC",steps=1)
select = stepwise(mod,direction="forward/backward",criterion="AIC",steps=2)
select = stepwise(mod,direction="forward/backward",criterion="AIC",steps=3)

select = stepwise(mod,direction="forward/backward",criterion="AIC")
select = stepwise(mod,direction="forward/backward",criterion="BIC")

## Accuracy for the best sub-models

observed = coffee$Localisation

### Accuracy values for best submodels
acc = rep(0,5) # Initialize a vector of accuracy values
for (k in 1:5) {
   select = stepwise(mod,direction="forward/backward",criterion="AIC",steps=k,trace=0)
   predictions = predict(select,type="class")
   acc[k] = mean(predictions==observed)
}   

### 10-fold cross-validated accuracy values for best submodels

cvacc = rep(0,5) # Initialize a vector of accuracy values

folds = fold(coffee,k=10,cat_col="Localisation")$".folds" # Create balanced segments
folds 

cvpredictions = rep("1",nrow(coffee)) # Initialize a vector of predicted classes

for (k in 1:5) {
   select = stepwise(mod,direction="forward/backward",criterion="AIC",steps=k,trace=0)
   for (j in 1:10) {
      train = coffee[folds!=j,]
      test = coffee[folds==j,]
      submod = multinom(formula(select),data=train,trace=FALSE,maxit=200) 
      cvpredictions[folds==j] = predict(submod,newdata=test,type="class")
   }
   cvacc[k] = mean(cvpredictions==coffee$Localisation)
}   

### Accuracy plot for stepwise feature selection
plot(1:5,acc,type="b",pch=16,xlab="Number of variables in the submodel",
     ylab="Accuracy",main="Quality of fit of best submodels",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,lwd=2,col="darkgray",bty="l")
lines(1:5,cvacc,type="b",col="blue",pch=15,lwd=2)
legend("topleft",lwd=2,pch=c(16,15),legend=c("Internal validation","Cross validation"),
       bty="n",cex=1.25,col=c("darkgray","blue"))
grid()

### 10-fold cross-validated accuracy for best submodel

folds = fold(coffee,k=10,cat_col="Localisation")$".folds" # Create balanced segments
folds 

cvpredictions = rep("1",nrow(coffee)) # Initialize a vector of predicted classes

for (j in 1:10) {
   train = coffee[folds!=j,]
   test = coffee[folds==j,]
   mod = multinom(Localisation~.,data=train,trace=FALSE,maxit=200) 
   select = stepwise(mod,direction="forward/backward",criterion="AIC",trace=0)
   cvpredictions[folds==j] = predict(select,newdata=test,type="class")
   print(paste("Segment ",j,sep=""))
}

mean(cvpredictions==coffee$Localisation)
