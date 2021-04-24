

# Set working directory

setwd("C:/Users/David/Dropbox/ADB_2021/Cours/Session1")

# Required packages

install.packages("nnet") # Installation is only needed if the package is missing

require(nnet)            # Multinomial logistic regression

# Regression with a real-valued response

## Import 'invasive probe' dataset

pig = read.table("./data/invasive.txt")
str(pig)  # Overview of data
dim(pig)  # Number of rows and columns in pig

summary(pig)   # Provides a columnwise summary of the data table

## Least-squares fit the most complete model

X = pig[,-1] # 60 x 11 matrix with measurements of explanatory variables   
y = pig$LMP  # 60-vector of response values

### Handmade calculations

Sx = var(X)               # Sample variance matrix
sxy = cov(X,y)            # Sample covariance
beta.calc = solve(Sx,sxy) # Estimated regression coefficients

### Compare with estimations using glm

mod = glm(LMP~.,data=pig) # LS fit of the model
data.frame(beta.calc=beta.calc,beta.mod=coef(mod)[-1])

## Assessment of the fit

### Using RSS and R2

RSS = sum(residuals(mod)^2)  # Residuals sum-of-squares

mod0 = glm(LMP~1,data=pig) # LS fit of the null model
RSS0 = sum(residuals(mod0)^2)  # Residuals sum-of-squares of null model
R2=(RSS0-RSS)/RSS0 # R2
R2

fitted.lmp = fitted(mod) # Fitted LMP values
observed.lmp = pig$LMP   # Observed LMP values
cor(observed.lmp,fitted.lmp)^2   # Same

### Scatterplot of fitted against observed LMP values

plot(observed.lmp,fitted.lmp,type="p",pch=16,bty="n",xlab="Observed LMP",
   ylab="Fitted LMP",main="Fitted versus observed LMP values",
   cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
abline(0,1)
text(50,65,paste("R2=",round(R2,3)),cex=1.25)

# Regression with a K-class response

## Import coffee data

coffee = read.table("./data/coffee.txt")
dim(coffee)       # Number of rows and columns in data
str(coffee)       # Overview of data
coffee$Localisation = factor(coffee$Localisation) # Convert 'Localisation' into a factor
summary(coffee)   # Provides a columnwise summary of the data table (8 first columns)

## ML fit the most complete model for 'Localisation'

mod = multinom(Localisation~.,data=coffee) # ML fit of the logistic model
coef(mod)

mod = multinom(Localisation~.,data=coffee,maxit=200) # ML fit of the logistic model
coef(mod)

### Assessment of the fit

deviance(mod)          # Residual deviance

mod0 = multinom(Localisation~1,data=coffee) # ML fit of the null model
deviance(mod0)-deviance(mod)                # Explained deviance

### Observed versus fitted values

proba = fitted(mod)            # Estimated probabilities of each class
observed_class = coffee$Localisation
head(data.frame(round(proba,3),observed_class))

# Confusion matrix

fitted_class = predict(mod,type="class")  # Bayes rule
head(fitted_class)

confusion = table(observed_class,fitted_class)
confusion

rowSums(confusion)
confusion_percentage = 100*confusion/outer(rowSums(confusion),rep(1,7))
round(confusion_percentage,3)
