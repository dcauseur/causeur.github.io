

# Set working directory

setwd("C:/Users/David/Dropbox/ADB_2021/Cours/Session3")

### Required packages

install.packages("glmnet") # Installation is only needed if the package is missing
install.packages("leaps")
install.packages("fields")
install.packages("pls")

require(leaps)             # For variable selection
require(glmnet)            # For penalized regression procedures
require(fields)            # For image.plot
require(pls)               # For cvsegments, plsr, ...

# Import 'scanner' dataset

pig = read.table("./Data/scanner.txt")
dim(pig)       # Numbers of rows and columns in dta
str(pig)       # Overview of the data table

# Data description

## Display curves

matplot(1:137,t(pig[,1:137]),type="l",lty=1,col="orange",
        main="CT curves for muscle percentage",xlab="Anatomical location",
        ylab="Muscle %")

## Display correlation curve

corxy = cor(pig[,1:137],pig[,138])
plot(1:137,corxy,type="l",col="orange",
        main="Correlation between CT curves for muscle percentage and LMP",
        xlab="Anatomical location",
        ylab="Correlation")

## Correlation matrix across explanatory variables
xcor = cor(pig[,-138])

### Image plot of the correlations across explanatory variables
image.plot(1:137,1:137,xcor[,137:1],xlab="Indices of explanatory variables",
           ylab="Indices of explanatory variables",
           main="Image plot of the correlations across explanatory variables",
           yaxt="n",cex.lab=1.25,cex.axis=1.25,cex.main=1.25)

# Unstability of least-squares fit on highly correlated variables

## Rank the models according to their RSS (using package leaps)
select = summary(regsubsets(LMP~.,data=pig,nvmax=115,method="forward"))
head(select$which)

## Coefficient plot for the best model with 5 variables

beta = rep(0,137)       # Initialize an empty vector of regression coefficients
mod = glm(LMP~.,data=pig[,c(select$which[5,-1],TRUE),drop=FALSE])
beta[select$which[5,-1]] = coef(mod)[-1]

plot(1:137,beta,bty="l",type="b",pch=16,lwd=2,xlab="Indices of explanatory variables",
     ylab="Regression coefficients",main="Regression coefficients",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
mtext("Best model with 5 variables",cex=1.25)

## Coefficient plot for the best model with 100 variables

beta = rep(0,137)      # Initialize an empty vector of regression coefficients
mod = glm(LMP~.,data=pig[,c(select$which[100,-1],TRUE),drop=FALSE])
beta[select$which[100,-1]] = coef(mod)[-1]

plot(1:137,beta,bty="l",type="b",pch=16,lwd=2,xlab="Indices of explanatory variables",
     ylab="Regression coefficients",main="Regression coefficients",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
mtext("Best model with 100 variables",cex=1.25)

# Ridge estimation of a regression model

## Matrix of explanatory variables
x = as.matrix(pig[,-138])

## Vector of response values
y = pig[,138]

## Decreasing sequence of log-lambda values
loglambda = seq(10,-10,length=100) 

## Fit ridge regression for a sequence of lambda
mod = glmnet(x,y,alpha=0,lambda=exp(loglambda)) 

## Initialize an empty matrix of predicted values
cvpred = matrix(0,nrow=nrow(pig),ncol=100) 

## Create 10 random segments of the dataset
segments = cvsegments(nrow(pig),10) 

for (k in 1:10) {
   mod = glmnet(x[-segments[[k]],],y[-segments[[k]]],alpha=0,lambda=exp(loglambda))
   cvpred[segments[[k]],] = predict(mod,newx=x[segments[[k]],]) 
}

## Predicted versus observed LMP with large lambda
plot(pig$LMP,cvpred[,1],pch=16,bty="n",xlab="Observed LMP",
     ylab="Predicted LMP",main="Fitted versus observed LMP",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,ylim=range(pig$LMP))
mtext(expression(Large~lambda),cex=1.25)
abline(0,1,lwd=2,col="darkgray")

## Predicted versus observed LMP with small lambda
plot(pig$LMP,cvpred[,100],pch=16,bty="n",xlab="Observed LMP",
     ylab="Predicted LMP",main="Fitted versus fitted LMP",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,ylim=range(pig$LMP))
mtext(expression(lambda~close~to~zero),cex=1.25)
abline(0,1,lwd=2,col="darkgray")

## Searching for the best penalty parameter
cvmod = cv.glmnet(x,y,alpha=0,lambda=exp(loglambda))

## MSEP profile in ridge regression
plot(cvmod)

## Take the predictions using the optimal lambda value
pred = cvpred[,which.min(cvmod$cvm)]

## Fitted versus observed LMP
plot(pig$LMP,pred,pch=16,bty="n",xlab="Observed LMP",
     ylab="Predicted LMP",main="Fitted versus fitted LMP",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,ylim=range(pig$LMP))
mtext(expression(Optimal~lambda),cex=1.25)
abline(0,1,lwd=2,col="darkgray")

## MSEP for the Ridge estimator

### Initialize an empty matrix of predicted values
cvpred = matrix(0,nrow=nrow(pig),ncol=100) 

### Create 10 random segments of the dataset
segments = cvsegments(nrow(pig),10) 

for (k in 1:10) {
   trainx = x[-segments[[k]],]
   trainy = y[-segments[[k]]]
   testx = x[segments[[k]],]
   cvmod = cv.glmnet(trainx,trainy,alpha=0,lambda=exp(loglambda))
   mod = glmnet(trainx,trainy,alpha=0,lambda=exp(loglambda))
   cvpred[segments[[k]],] = predict(mod,newx=testx)[,which.min(cvmod$cvm)] 
   print(paste("Segment ",k,sep=""))
}

MSEP = mean((pig$LMP-cvpred)^2)

## Searching for the best penalty parameter
cvmod = cv.glmnet(x,y,alpha=0,lambda=exp(loglambda))

## MSEP profile in ridge regression
plot(cvmod,ylim=c(1,5))
abline(h=MSEP,col="orange",lwd=2)

# LASSO estimation of a regression model

## Univariate model

x = pig$X92
y = pig$LMP

n = nrow(pig) ; n
sxy = cov(x,y) ; sxy
s2x = var(x) 

vlambda = seq(0,150,length=1000)
beta = (sxy-vlambda/(2*n))/s2x 
beta[sxy<=vlambda/(2*n)] = 0

# Estimated regression coefficient along with lambda
plot(vlambda,beta,type="l",lwd=2,bty="n",xlab=expression(lambda),
     ylab=expression(hat(beta)),main="Estimated regression coefficient using LASSO",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,col="darkgray")

## Multivariate model

### Matrix of explanatory variables
x = as.matrix(pig[,-138])

### Vector of response values
y = pig[,138]

### Decreasing sequence of log-lambda values
loglambda = seq(10,-10,length=100) 

### Fit LASSO regression for a sequence of lambda
mod = glmnet(x,y,alpha=1,lambda=exp(loglambda)) 

### Initialize an empty matrix of predicted values
cvpred = matrix(0,nrow=nrow(pig),ncol=100) 

### Create 10 random segments of the dataset
segments = cvsegments(nrow(pig),10) 

for (k in 1:10) {
   mod = glmnet(x[-segments[[k]],],y[-segments[[k]]],alpha=1,lambda=exp(loglambda))
   cvpred[segments[[k]],] = predict(mod,newx=x[segments[[k]],]) 
}

cvmod = cv.glmnet(x,y,alpha=1,lambda=exp(loglambda))
### MSEP profile in ridge regression
plot(cvmod)

### Take the predictions using the optimal lambda value
pred = cvpred[,which.min(cvmod$cvm)]

### Fitted versus observed LMP
plot(pig$LMP,pred,pch=16,bty="n",xlab="Observed LMP",
     ylab="Predicted LMP",main="Fitted versus fitted LMP",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25,ylim=range(pig$LMP))
mtext(expression(Optimal~lambda),cex=1.25)
abline(0,1,lwd=2,col="darkgray")

### Sparsity of the estimated regression model
mod = glmnet(x,y,alpha=1,lambda=exp(loglambda))
plot(1:137,mod$beta[,which.min(cvmod$cvm)],pch=16,bty="n",xlab="Indices of explanatory variables",
     ylab=expression(hat(beta)),main="Estimated regression coefficients using LASSO",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
mtext(expression(Optimal~lambda),cex=1.25)

### MSEP for the Lasso estimator

#### Initialize an empty matrix of predicted values
cvpred = matrix(0,nrow=nrow(pig),ncol=100) 

#### Create 10 random segments of the dataset
segments = cvsegments(nrow(pig),10) 

for (k in 1:10) {
   trainx = x[-segments[[k]],]
   trainy = y[-segments[[k]]]
   testx = x[segments[[k]],]
   cvmod = cv.glmnet(trainx,trainy,alpha=1,lambda=exp(loglambda))
   mod = glmnet(trainx,trainy,alpha=1,lambda=exp(loglambda))
   cvpred[segments[[k]],] = predict(mod,newx=testx)[,which.min(cvmod$cvm)] 
   print(paste("Segment ",k,sep=""))
}

MSEP = mean((pig$LMP-cvpred)^2)

## Searching for the best penalty parameter
cvmod = cv.glmnet(x,y,alpha=1,lambda=exp(loglambda))

## MSEP profile in lasso regression
plot(cvmod,ylim=c(1,5))
abline(h=MSEP,col="orange",lwd=2)

# Lasso estimation of a multinomial logistic regression model

## Import data

coffee_nirs = read.table("./Data/coffee_nirs.txt")
coffee_nirs$Localisation = factor(coffee_nirs$Localisation)
dim(coffee_nirs)
str(coffee_nirs[,1:15])

## Standard Normal Variate transformation of NIRS

x = coffee_nirs[,-(1:6)]    # NIRS data
snv_x = t(scale(t(x)))      # SNV-transformed NIRS 

## Displays the NIRS

wn = seq(402,2500,2)           # Wavenumbers
matplot(wn,t(snv_x),type="l",lwd=2,col="orange",lty=1,
        bty="l",xlab=expression(Wave~numbers~(nm^-1)),
        ylab="SNV-transformed NIRS",main="NIRS data of coffee samples",
        cex.main=1.25,cex.lab=1.25,cex.axis=1.25)

## LASSO estimation of the regression model

y = coffee_nirs$Localisation

### Choice of the best penalty parameter

loglambda = seq(-10,-1,length=100)

coffee_nirs.cvlasso = cv.glmnet(snv_x,y,family="multinomial",type.measure="deviance",
                                lambda=exp(loglambda))

### CV'd residual deviance for each penalty parameter

plot(coffee_nirs.cvlasso)

### Assessment of the model obtained with optimal lambda

coffee_nirs.lasso = glmnet(snv_x,y,family="multinomial",lambda=exp(loglambda))

proba = predict(coffee_nirs.lasso,newx=snv_x,
                type="response")[,,which.min(coffee_nirs.cvlasso$cvm)]
head(proba)   # For each coffee, estimated class probabilities

predictions = predict(coffee_nirs.lasso,newx=snv_x,
                type="class")[,which.min(coffee_nirs.cvlasso$cvm)]
head(predictions)   # For each coffee, predicted class using Bayes rule

confusion = table(coffee_nirs$Localisation,predictions)

acc = mean(coffee_nirs$Localisation==predictions)
acc

### Accuracy of the Lasso regression model

dta = data.frame(snv_x,"Localisation"=y)
dim(dta)

segs = fold(dta,k=10,cat_col="Localisation")$".folds"
# 10-fold balanced partition of the sample

cvpredictions = rep("0",nrow=nrow(dta))
# Empty matrix to store the CV'd probabilities,
# for each coffee, of being located in Loc. 6

for (k in 1:10) {
   train = dta[segs!=k,]
   test = dta[segs==k,]
   dta.cvlasso = cv.glmnet(as.matrix(train[,-1051]),train[,1051],family="multinomial",
                           type.measure="deviance",lambda=exp(loglambda))
   dta.lasso = glmnet(as.matrix(train[,-1051]),train[,1051],family="multinomial",
                           lambda=exp(loglambda))
   cvpredictions[segs==k] = predict(dta.lasso,newx=as.matrix(test[,-1051]),
                              type="class")[,which.min(dta.cvlasso$cvm)]
   print(paste("Segment ",k," over 10",sep=""))
}

cv.acc = mean(coffee_nirs$Localisation==cvpredictions)
cv.acc # 0.95







