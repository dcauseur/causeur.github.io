

# Set working directory

setwd("C:/Users/David/Dropbox/ADB_2021/Cours/Session4")

### Required packages

install.packages("glmnet") # Installation is only needed if the package is missing
install.packages("leaps")
install.packages("fields")
install.packages("pls")
install.packages("MASS")
install.packages("viridis")
install.packages("groupdata2")

require(leaps)             # For variable selection
require(glmnet)            # For penalized regression procedures
require(fields)            # For image.plot
require(pls)               # For cvsegments, plsr, ...
require(MASS)              # For LDA
require(viridis)           # For nice colors
require(groupdata2)        # For balanced cross-validation

# Import 'scanner' dataset

pig = read.table("./Data/scanner.txt")
dim(pig)       # Numbers of rows and columns in dta
str(pig)       # Overview of the data table

# PLS with one latent component

## Matrix of explanatory variables
x = pig[,-138]

## Matrix of scaled explanatory variables
xstar = scale(x)

## Vector of response values
y = pig[,138]

## Squared covariances between scaled xs and y
s2xy = cov(xstar,y)^2
max(s2xy) 

## PLS fit with one component
lmp.pls = plsr(LMP~.,data=pig,ncomp=1,scale=TRUE)

## Extract latent variable and corresponding loadings
lv = scores(lmp.pls)
alpha = loadings(lmp.pls)

## Squared covariance between scaled latent variable and y 
s2ty = cov(scale(lv[,1]),y)^2
s2ty

## loadings plot for the PLS model with one component
plot(1:137,alpha[,1],bty="l",type="b",pch=16,lwd=2,xlab="Indices of explanatory variables",
     ylab="Loadings",main="Loadings of the best latent variable",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)

## Latent variable modeling

mod = lm(LMP~LV,data=data.frame(LMP=y,LV=lv[,1]))
fitted(mod)[1:6]

fitted(lmp.pls)[1:6,1,1]

# RMSEP of the model
lmp.pls = plsr(LMP~.,data=pig,ncomp=1,scale=TRUE,validation="CV",segments=10)
RMSEP(lmp.pls)

# Percentage of explained variance
explvar(lmp.pls)

### Adding a 2nd latent variable
lmp.pls2 = plsr(LMP~.,data=pig,ncomp=2,scale=TRUE)
lv2 = scores(lmp.pls2)
cor(lv2)

# Percentage of explained variance
explvar(lmp.pls2)

# Score plot
scoreplot(lmp.pls2,pch=16,ylim=range(lv2[,1]))

# Does it improve the fit?
R2(lmp.pls2)

# Does it improve the prediction accuracy?
lmp.pls2 = plsr(LMP~.,data=pig,ncomp=2,scale=TRUE,validation="CV",segments=10)
RMSEP(lmp.pls2)

### How many latent variables to introduce in the model?
lmp.pls = plsr(LMP~.,data=pig,ncomp=100,scale=TRUE,validation="CV",segments=10)

# Percentage of explained variance
explvar(lmp.pls)

# Selecting the optimal number of latent variables
selectNcomp(lmp.pls,method="onesigma",plot=TRUE)

### Implementation of a complete 10-fold CV procedure 

n = nrow(pig)
cvpred = rep(0,n)

segs = cvsegments(n,10)
cvpred = rep(0,n)

for (k in 1:10) {
   cvmod = plsr(LMP~.,data=pig[-segs[[k]],],ncomp=30,scale=TRUE,validation="CV",segments=10)
   bestncomp = selectNcomp(lmp.pls,method="onesigma")
   cvpred[segs[[k]]] = predict(cvmod,newdata=pig[segs[[k]],])[,,bestncomp]
   print(k)
}

PRESS.pls = sum((pig$LMP-cvpred)^2) ; PRESS.pls

# Linear Discriminant Analysis 

# Import oyster data 

oysters = read.table("./Data/oysters.txt",stringsAsFactors=TRUE)
str(oysters)

# Standard notations
n = nrow(oysters)
p = ncol(oysters)-1
K = length(table(oysters$Sp))

# Prior probabilities
priors = table(oysters$Sp)/n
round(priors,3)

# Class means
mu = tapply(oysters$C2,INDEX=oysters$Sp,FUN=mean)

# Within class standard deviation
sigma = summary(lm(C2~Sp,data=oysters))$sigma

## Plot of posterior probabilities

vec_C2 = seq(20,60,length=1000) # A sequence of C2 values

f = outer(vec_C2,mu,dnorm,sd=sigma)   # All values of f(x,mu_k,sigma)
denom = (f%*%priors)[,1]            
posteriors = f*outer(1/denom,priors,"*")

matplot(vec_C2,posteriors,type="l",bty="l",lty=1,lwd=2,
        xlab="C2",ylab="Posterior class probabilities",
        main="Posterior probabilities of being of a given species",
        cex.lab=1.25,cex.axis=1.25,cex.main=1.25,col=viridis(7))
legend(40,1,col=viridis(7),lwd=2,lty=1,legend=levels(oysters$Sp),
       bty="n",cex=1.25)

## Prediction on the learning dataset

f = outer(oysters$C2,mu,dnorm,sd=sigma)   # All values of f(x,mu_k,sigma)
denom = (f%*%priors)[,1]            
posteriors = f*outer(1/denom,priors,"*")

whichmap = apply(posteriors,1,which.max)
prediction = levels(oysters$Sp)[whichmap]
# Which class has maximal posterior probability?

table(oysters$Sp,prediction)

## Fisher's LDA score using lda

oysters.lda = lda(Sp~.,data=oysters)
beta = coef(oysters.lda)[,1]    # Fisher's linear coefficients
L = predict(oysters.lda)$x[,1]  # Fisher's score

## Separation between classes

### F-test statistic for the linear score
anova(lm(L~oysters$Sp)) # Note: within-class variance of L is 1

oysters.lda$svd[1]^2

### F-test statistics for each explanatory variable
Ftest = rep(0,p)
for (k in 1:p) {
   Ftest[k] = anova(lm(C~Sp,data=data.frame(C=oysters[,k],Sp=oysters$Sp)))[1,4]
}
Ftest
   
### Prediction based on linear score - Geometrical interpretation

predictions.map = predict(oysters.lda,dimen=1)$class

L.classmeans = tapply(L,INDEX=oysters$Sp,FUN=mean)
class.distances = outer(L,L.classmeans,"-")
head(class.distances)

predictions.mindist = apply(abs(class.distances),1,which.min)
predictions.mindist = levels(oysters$Sp)[predictions.mindist]
   
table(predictions.map,predictions.mindist)

## Complete K-class LDA

oysters.lda$svd^2   # F-statistics
scores = predict(oysters.lda)$x

colors = as.numeric(oysters$Sp)
plot(oysters.lda,dimen=2,col=colors,bty="n",cex=1.25,cex.axis=1.25,
     cex.lab=1.25,cex.main=1.25,main="First two LD scores")
abline(h=0,v=0,col="darkgray",lwd=2)

## Number of LD scores in prediction

### Prediction on the training sample

predictions.map = predict(oysters.lda,dimen=2)$class
table(oysters$Sp,predictions.map,dnn=list("observed","predicted"))
mean(oysters$Sp==predictions.map) # accuracy

predictions.map = predict(oysters.lda,dimen=6)$class
table(oysters$Sp,predictions.map,dnn=list("observed","predicted"))
mean(oysters$Sp==predictions.map) # accuracy

### Leave-one-out prediction of the MAP (6 linear scores) 

oysters.lda = lda(Sp~.,data=oysters,CV=TRUE)
predictions.map = oysters.lda$class
table(oysters$Sp,predictions.map,dnn=list("observed","predicted"))
mean(oysters$Sp==predictions.map) # Leave-one-out accuracy

