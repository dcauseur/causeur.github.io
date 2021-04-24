

# Set working directory

setwd("C:/Users/David/Dropbox/ADB_2021/Cours/Session5")

# Required packages

install.packages("gam") # Installation is only needed if the package is missing
install.packages("pls")
install.packages("fields")

require(gam)             # For generalized additive models
require(pls)             # For cvsegments ...
require(fields)          # For image.plot

# Import 'ozone' dataset

ozone = read.table("./data/ozone.txt",sep=";",header=TRUE,row.names=1)
str(ozone) # Overview of Ozone data 

# Ozone concentration along temperature

## Scatterplot for maximal ozone concentration along temperature at 12PM
plot(maxO3~T12,data=ozone,bty="l",type="p",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)

## k-neighborhood of T12=28
n = nrow(ozone)
k = round(0.20*n)          # 20% of the sample in the neighborhood
dx0 = abs(ozone$T12-28)
dx0.sorted = sort(dx0)     # Sort in ascending order
Nx0 = (1:n)[which(dx0<=dx0.sorted[k])] # Indices of individuals within the k-neighborhood of T12=28
lim.inf = min(ozone$T12[Nx0])
lim.inf                                # Smallest T12 value in the k-neighborhood of T12=28
lim.sup = max(ozone$T12[Nx0])
lim.sup                                # Largest T12 value in the k-neighborhood of T12=28

## k-neighborhood of T12=28
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
polygon(c(lim.inf,lim.sup,lim.sup,lim.inf),c(30,30,170,170),col="papayawhip")
points(ozone$T12,ozone$maxO3,type="p",pch=16)
abline(v=28,lwd=2,col="orange")

## Tricube weights within N(x0=28)
u = rep(1,n)           # Initialize a vector of u values with 1 for all
u[Nx0] = dx0[Nx0]/max(dx0[Nx0])       # Within the neighborhood, the values of u are changed
omega = (1-u^3)^3                     # Tricube values
ord = order(ozone$T12)                # Gives the indices of individuals for an ascending sorting of T12

## Display the tricube weights within N(x0=28)
plot(ozone$T12[ord],omega[ord],bty="l",type="n",lwd=2,xlab="Temperature at 12h",
     ylab="Weights",main="Weights for local polynomial approximation at T12=28",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
polygon(c(lim.inf,lim.sup,lim.sup,lim.inf),c(0,0,1,1),col="papayawhip")
lines(ozone$T12[ord],omega[ord],type="l",lwd=2)
abline(v=28,lwd=2,col="orange")

## Weighted local polynomial approximation in N(x0=28)
local.poly = lm(maxO3~T12+I(T12^2),data=ozone,weights=omega)
prediction = predict(local.poly,newdata=data.frame(T12=28))
prediction

seq.t12 = seq(lim.inf,lim.sup,length=1000)
fit = predict(local.poly,newdata=data.frame(T12=seq.t12))

## Local polynomial approximation in N(T12=28)
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
polygon(c(lim.inf,lim.sup,lim.sup,lim.inf),c(30,30,170,170),col="papayawhip")
points(ozone$T12,ozone$maxO3,type="p",pch=16)
abline(v=28,lwd=2,col="orange")
lines(seq.t12,fit,lwd=5,col="darkgray")

# Generalization to a sequence of T12 values: loess fit using gam
t12.loess = gam(maxO3~lo(T12,span=0.2),data=ozone)        # Gam fit with 20%-neighborhoods
seq.t12 = seq(5,35,length=1000)                           # A grid of T12 values 
fit = predict(t12.loess,newdata=data.frame(T12=seq.t12))  # Fitted values

# Local polynomial approximation with 20%-neighborhoods
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

# Local polynomial approximation with 5%-neighborhoods
t12.loess = gam(maxO3~lo(T12,span=0.95),data=ozone)
seq.t12 = seq(5,35,length=1000)
fit = predict(t12.loess,newdata=data.frame(T12=seq.t12))

plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

## Choice of k using cross-validation
train = sample(1:nrow(ozone),800)   # Indices of individuals in the training dataset
test = setdiff(1:nrow(ozone),train) # Indices of individuals in the test dataset
ozone.train = ozone[train,]
ozone.test = ozone[test,]

seq.span = seq(0.05,0.95,length=100)  # A grid of span values
PRESS = rep(0,length(seq.span))       # Initialize a vector of PRESS statistics 

for (i in 1:length(seq.span)) {       # Cycle over span values
   t12.loess = gam(maxO3~lo(T12,span=seq.span[i]),data=ozone.train)
   pred = predict(t12.loess,newdata=ozone.test)
   PRESS[i] = sum((ozone.test$maxO3-pred)^2)
   print(paste(i,"th span value",sep=""))
}

## PRESS plot to choose k
plot(seq.span,PRESS,bty="l",type="l",lwd=2,xlab="Bandwidth in loess",
     ylab="PRESS",main="Choice of the bandwidth in loess using CV",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
abline(v=seq.span[which.min(PRESS)],lwd=2,col="orange")

## Fit GAM with optimal span
t12.loess = gam(maxO3~lo(T12,span=seq.span[which.min(PRESS)]),data=ozone)
seq.t12 = seq(5,35,length=1000)
fit = predict(t12.loess,newdata=data.frame(T12=seq.t12))

## Local polynomial approximation with optimal k
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

# Spline smoothing

## B-spline basis (degree=1)
knots = seq(5,35,length=7)
basis = bs(ozone$T12,knots=knots[-c(1,length(knots))],degree=1,
           intercept=TRUE)
ord = order(ozone$T12)

## Display B-splines of degree 1
matplot(ozone$T12[ord],basis[ord,],bty="l",type="l",lwd=2,lty=1,
        xlab="Temperature at 12h",
        ylab="B-spline basis",main="B-spline basis of degree 1",cex.lab=1.25,
        cex.axis=1.25,cex.main=1.25)

## Cubic B-spline basis (degree=3)
basis = bs(ozone$T12,knots=knots[-c(1,length(knots))],intercept=TRUE)

## Display Cubic B-splines
matplot(ozone$T12[ord],basis[ord,],bty="l",type="l",lwd=2,lty=1,
        xlab="Temperature at 12h",
        ylab="B-spline basis",main="Cubic B-spline basis",cex.lab=1.25,
        cex.axis=1.25,cex.main=1.25)

## Spline smoothing with 5 interior knots
t12.spline = lm(maxO3~-1+.,data=data.frame(maxO3=ozone$maxO3,basis))
coef(t12.spline)[1:6]

seq.t12 = seq(5,35,length=1000)  # A grid of T12 values
basis = bs(seq.t12,knots=knots[-c(1,length(knots))],intercept=TRUE) # Corresponding B-spline basis 
fit = predict(t12.spline,newdata=data.frame(basis))

## Display spline approximation (with 5 interior knots)
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

## Spline smoothing with 18 interior knots
knots = seq(5,35,length=20)
basis = bs(ozone$T12,knots=knots[-c(1,length(knots))],intercept=TRUE)
t12.spline = lm(maxO3~-1+.,data=data.frame(maxO3=ozone$maxO3,basis))
basis = bs(seq.t12,knots=knots[-c(1,length(knots))],intercept=TRUE)
fit = predict(t12.spline,newdata=data.frame(basis))

## Display spline approximation (with 18 interior knots)
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

## Nonparametric degree of freedom
t12.spline = gam(maxO3~s(T12,df=5),data=ozone) # Arbitrarily, 5 degrees of freedom
seq.t12 = seq(5,35,length=1000)                # A grid of T12 values
fit = predict(t12.spline,newdata=data.frame(T12=seq.t12))

## Spline approximation
plot(ozone$T12,ozone$maxO3,bty="l",type="n",pch=16,xlab="Temperature at 12h",
     ylab="Maximal ozone concentration",
     main="Maximal ozone concentration along the temperature at 12h",cex.lab=1.25,
     cex.axis=1.25,cex.main=1.25)
points(ozone$T12,ozone$maxO3,type="p",pch=16)
lines(seq.t12,fit,lwd=5,col="orange")

t12.loess = gam(maxO3~lo(T12,span=0.35),data=ozone)
summary(t12.loess)$df
fit = predict(t12.loess,newdata=data.frame(T12=seq.t12))

## Adds the loess approximation with similar number of degrees of freedom
lines(seq.t12,fit,lwd=5,col="blue")

# Nonparametric ANOVA

t12.loess = gam(maxO3~lo(T12,span=0.35),data=ozone)
t12.lm = gam(maxO3~T12,data=ozone)

anova(t12.lm,t12.loess,test="F")

anova(t12.loess)

# Prediction accuracy

n = nrow(ozone)              # Sample size

segments = cvsegments(n,10)  # 10 random segments
cvpred.lo = rep(0,n)         # Initialize a n-vector of CV'd predictions using loess
cvpred.lm = rep(0,n)         # Initialize a n-vector of CV'd predictions using lm

for (k in 1:10) {
   dta.lo = gam(maxO3~lo(T12,span=0.35),data=ozone[-segments[[k]],])
   dta.lm = gam(maxO3~T12,data=ozone[-segments[[k]],])
   cvpred.lo[segments[[k]]] = predict(dta.lo,newdata=ozone[segments[[k]],])
   cvpred.lm[segments[[k]]] = predict(dta.lm,newdata=ozone[segments[[k]],])
   print(paste(k,"th segment",sep=""))
}

PRESS.lo = sum((ozone$maxO3-cvpred.lo)^2) ; PRESS.lo
PRESS.lm = sum((ozone$maxO3-cvpred.lm)^2) ; PRESS.lm

100*(PRESS.lm-PRESS.lo)/PRESS.lm   # Gain percentage by using loess 

# Multivariate additive models

## Model fitting

dta.gam = gam(maxO3~lo(T6,span=0.35)+lo(T12,span=0.35)+lo(T15,span=0.35)+lo(Ne6,span=0.35)+
                 lo(Ne12,span=0.35)+lo(Ne15,span=0.35)+lo(Vx,span=0.35)+lo(maxO3v,span=0.35),
              data=ozone)

# Plot of marginal effects
par(mfrow=c(2,4))
plot(dta.gam,lwd=3,col="blue",cex.lab=1.25,se=TRUE)
par(mfrow=c(1,1))

## Nonparametric anova

anova(dta.gam)

## Prediction performance of the full gam 

segments = cvsegments(n,10)   # Random 10-segments
cvpred.gam = rep(0,n)         # Initialize a n-vector of CV'd predicted values using gam

for (k in 1:10) {
   dta.gam = gam(maxO3~lo(T6,span=0.35)+lo(T12,span=0.35)+lo(T15,span=0.35)+lo(Ne6,span=0.35)+
                    lo(Ne12,span=0.35)+lo(Ne15,span=0.35)+lo(Vx,span=0.35)+
                    lo(maxO3v,span=0.35),data=ozone[-segments[[k]],])
   cvpred.gam[segments[[k]]] = predict(dta.gam,newdata=ozone[segments[[k]],])
   print(paste(k,"th segment",sep=""))
}

PRESS.gam = sum((ozone$maxO3-cvpred.gam)^2) ; PRESS.gam

100*(PRESS.lo-PRESS.gam)/PRESS.lo # Gain in prediction performance w.r.t model with T12 only

### Stepwise model selection

O3.gam = gam(maxO3~1,data=ozone)
O3.step = step.Gam(O3.gam,scope=list("T15"=~1+T15+lo(T15,0.35),
                                     "maxO3v"=~1+maxO3v+lo(maxO3v,0.35),
                                     "T6"=~1+T6+lo(T6,0.35),
                                     "T12"=~1+T12+lo(T12,0.35),
                                     "Ne6"=~1+Ne6+lo(Ne6,0.35),
                                     "Ne12"=~1+Ne12+lo(Ne12,0.35),
                                     "Ne15"=~1+Ne15+lo(Ne15,0.35),
                                     "Vx"=~1+Vx+lo(Vx,0.35)))

segments = cvsegments(n,10)  # Random 10-segments of data
cvpred.stepgam = rep(0,n)    # Initialize a n-vector of CV'd predictions

for (k in 1:10) {
   dta.gam = gam(formula(O3.step),data=ozone[-segments[[k]],])
   cvpred.stepgam[segments[[k]]] = predict(dta.gam,newdata=ozone[segments[[k]],])
   print(paste(k,"th segment",sep=""))
}

PRESS.stepgam = sum((ozone$maxO3-cvpred.stepgam)^2) ; PRESS.stepgam

100*(PRESS.gam-PRESS.stepgam)/PRESS.gam # Gain in prediction performance w.r.t full model
