
### Set working directory

setwd("c:/users/david/dropbox/sad2020/session3/data")

### Required packages

install.packages("nnet") # Only if needed
require(nnet)            # For multinomial logistic regression

install.packages("questionr") # Only if needed
require(questionr)            # For odds-ratio

install.packages("aod") # Only if needed
require(aod)            # For splitbin

install.packages("Deriv") # Only if needed
require(Deriv)            # For differentiation

### Import data

# Import data stored in file maturity.txt
dta = read.table("maturity.txt")
# Convert the Maturity variable into a 3-level factor
dta$Maturity = factor(dta$Maturity,levels=c('1','2','3'))
# Overview of data
str(dta)

### 2-group maturity

dtab = subset(dta, Variety == "B" & Maturity != "3")
# Keeps only rows with both Variety == "B" and Maturity != "3"
dtab = droplevels(dtab)
# Forgets the unused levels

### Likelihood

# Let us take the first apricot on the sample
dtab[1,]

# Suppose beta0=-5, beta1=-1, to which extent is it likely?
beta0 = -5
beta1 = 0.2

# Linear score for the apricot
score = beta0+beta1*dtab$a[1]
score

# Probability that the apricot is in Maturity group 2 (Y=+1)
1/(1+exp(-score))

# Same with the last apricot in the sample
dtab[79,]

# Linear score for the apricot
score = beta0+beta1*dtab$a[79]
score

# Probability that the apricot is in Maturity group 1 (Y=-1)
1-(1/(1+exp(-score)))
1/(1+exp(score))

# Now taking all the apricots together

beta0 = -5
beta1 = 0.2
y = ifelse(dtab$Maturity=="2",1,-1)
scores = beta0+beta1*dtab$a
probas = 1/(1+exp(-y*scores))
likelihood = prod(probas)
likelihood

# Same with deviance

-2*log(likelihood)

### Illustration of Fisher scoring algorithm

## Let us suppose that beta0 = -5

# Deviance function
deviance_reglog = function(beta1,a,y,beta0=-5) {
   scores = beta0+beta1*a
   2*sum(log(1+exp(-y*scores)))
}

# Illustrative example 
deviance_reglog(beta1=0.2,a=dtab$a,y=y)

# Plot of the deviance function
b1 = seq(from=0,to=1,length=1000)
dev = rep(0,1000)
for (k in 1:1000) dev[k] = deviance_reglog(beta1=b1[k],a=dtab$a,y=y)

plot(b1,dev,xlab=expression(beta[1]),ylab="Deviance",lwd=2,col="orange",
     type="l",main="Deviance function",cex.lab=1.25,cex.axis=1.25)

# 1st order derivative of the deviance function
diff_deviance = function(beta1,a,y,beta0=-5) {
   scores = beta0+beta1*a
   probas = 1/(1+exp(-y*scores))
   resids = y*(1-probas)
   -2*sum(a*resids)
}

# Plot of the first order derivative
diffdev = rep(0,1000)
for (k in 1:1000) diffdev[k] = diff_deviance(beta1=b1[k],a=dtab$a,y=y)

plot(b1,diffdev,xlab=expression(beta[1]),ylab="Deviance",lwd=2,col="orange",
     type="l",main="Deviance function",cex.lab=1.25,cex.axis=1.25)
abline(h=0)

## Newton-Raphson algorithm

# 2nd order derivative of the deviance function
diff2_deviance = function(beta1,a,y,beta0=-5) {
        scores = beta0+beta1*a
        pi = 1/(1+exp(-scores))
        ve = pi*(1-pi)
        2*sum(a^2*ve)
}

# Let us start from beta1=0.3
# Plot of the tangent line at beta1
x0 = 0.3
y0 = diff_deviance(beta1=0.3,a=dtab$a,y=y)
slope = diff2_deviance(beta1=0.3,a=dtab$a,y=y)
abline(a=y0-slope*x0,b=slope,lty=5,col="darkgray",lwd=2)

# beta1 is updated
beta1 = x0-y0/slope
points(beta1,0,pch=16,cex=2,col="blue")

# Plot of the tangent line at beta 1
x0 = beta1
y0 = diff_deviance(beta1=beta1,a=dtab$a,y=y)
slope = diff2_deviance(beta1=beta1,a=dtab$a,y=y)
abline(a=y0-slope*x0,b=slope,lty=5,col="darkgray",lwd=2)

# beta1 is updated
beta1 = x0-y0/slope
points(beta1,0,pch=16,cex=2,col="blue")

# Now the whole iterative algorithm (starting with beta1=0)
maturity.logit = glm(Maturity~a,data=dtab,
                     family=binomial(link=logit),trace=TRUE,maxit=25,epsilon=1e-08)
# trace=TRUE prints the deviance along the iterations   

# Extract estimated regression parameters
beta = coef(maturity.logit)
beta

# Extract minimal deviance (residual deviance)
deviance(maturity.logit)

# Wald's z-tests for the significance of parameters
summary(maturity.logit)$coefficients

# Confidence intervals with confidence level 0.95
confint.default(maturity.logit)

