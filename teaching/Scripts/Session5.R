
setwd("C:/Users/David/Dropbox/DS20/Session5")

# Import data

dta = read.table("pig.txt",header=TRUE,stringsAsFactors=TRUE)
str(dta)  # Provides a columnwise overview of the data table

# Least-squares fit of the regression model
lmp.lm = lm(LMP~BFAT,data=dta)

# Confidence intervals of coefficients

## Confidence intervals for the regression coefficients 
## level = 0.95 (default) sets the confidence level at 0.95
cbind(coef(lmp.lm),confint(lmp.lm,level=0.95)) 

# Confidence band for regression line

x = c(12,23) # Two arbitrary values for fat depths
## Calculates corresponding fitted value
## with 95%-confidence intervals
pred = predict(lmp.lm,newdata=data.frame(BFAT=x),interval="confidence")
## pred has 3 columns: fit, upper and lower limits 
pred

## Calculates corresponding predictions
## with 95%-confidence intervals
pred = predict(lmp.lm,newdata=data.frame(BFAT=x),interval="prediction")
pred

# Assessment of the model fit

## Graphical assessment
fit = fitted(lmp.lm) # Extracts fitted values
plot(dta$LMP,fit,pch=16,xlim=c(48,67),ylim=c(48,67),
     xlab="Observed LMP values",ylab="Fitted LMP values")
abline(a=0,b=1,lwd=2,col="gray") # Adds the line y=x to the plot

## R2
summary(lmp.lm)$r.squared
cor(dta$LMP,fit)^2
cor(dta$LMP,dta$BFAT)^2

## F-test statistic
summary(lmp.lm)$fstatistic

# ANOVA table

anova(lmp.lm)

# Comparing regression lines

lmp.lm = lm(LMP~BFAT*GENET,data=dta)
coef(lmp.lm)

lmp.lm0 = lm(LMP~BFAT+GENET,data=dta)
coef(lmp.lm0)

RSS = sum(residuals(lmp.lm)^2) # Sum of squared residuals of the full model
RSS0 = sum(residuals(lmp.lm0)^2) # Sum of squared residuals of the submodel

Fstat = ((RSS0-RSS)/2)/(RSS/(59-6))
Fstat

anova(lmp.lm0,lmp.lm)
