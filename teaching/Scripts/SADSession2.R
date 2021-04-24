
### Set working directory

setwd("c:/users/david/dropbox/sad2020/session2/data")

### Required packages

install.packages("nnet") # Only if needed
require(nnet)            # For multinomial logistic regression

install.packages("questionr") # Only if needed
require(questionr)            # For odds-ratio

install.packages("aod") # Only if needed
require(aod)            # For splitbin

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
summary(dtab$a)
# Summary statistics for 'a'
table(dtab$Maturity)
# Counts apricots in each maturity group

# Partition of the 'a' values
vbreaks = seq(from=-8,to=24,by=2) 
# Counts in each bin of the partition 
hista = hist(x=dtab$a,breaks=vbreaks,plot=FALSE)
hista$counts

# Counts only the apricots in maturity level 2
hista2 = hist(x=dtab$a[dtab$Maturity=="2"],breaks=vbreaks,plot=FALSE)
hista2$counts

plot(hista$mids,hista2$counts/hista$counts,type="b",pch=16,xlab="a values",
     ylab="Proportion of apricots in maturity level 2",lwd=2,
     main="Empirical effect curve of a on the maturity")
mtext("Variety B")

# Fits the logit regression model
maturity.logit = glm(Maturity~a,data=dtab,family=binomial(link=logit))
# Extracts estimated coefficients
beta = coef(maturity.logit)
beta

# veca is a high resolution sequence of 'a' values
veca = seq(from=min(dtab$a),to=max(dtab$a),length=1000)

# Linear score 
linear_score = beta[1]+beta[2]*veca
head(linear_score)

# Linear score (same with function predict)
linear_score = predict(maturity.logit,newdata=data.frame(a=veca))
head(linear_score)

# Probability of being in maturity level 2 
proba = 1/(1+exp(-linear_score)) # Inverse-logit transformation
head(proba)

# Probability of being in maturity level 2 (same with predict)
proba = predict(maturity.logit,newdata=data.frame(a=veca),type="response")
head(proba)

# Adds to the plot
lines(veca,proba,lwd=2,col="orange")

## Odds-ratio of color-index a

odds.ratio(maturity.logit)

### 3-group maturity

# Fits the logit regression model
maturity.mlogit = multinom(Maturity~a,data=dta)
# Extracts estimated coefficients
beta = coef(maturity.mlogit)
beta

proba = predict(maturity.mlogit,newdata=data.frame(a=veca),type="probs")
head(proba)

# Adds to the plot
matplot(veca,proba,lwd=2,type="l",col=c("coral","coral2","coral4"),xlab="a values",
        main="Probability curves for 3-group maturity",ylab="Probability",lty=1,
        ylim=c(0,1))
legend("topright",lwd=2,lty=1,bty="n",col=c("coral","coral2","coral4"),
       legend=c("Group 1","Group 2","Group 3"))

### Two-way analysis of variance for the lamb study

# Reporting the experimental results into the data frame 'lamb'
lamb = data.frame(Feeding=c("F1","F1","F2","F2"),Housing=c("H1","H2","H1","H2"),
                  Coloured=c(10,12,15,16),Total=rep(20,times=4))
lamb

# From aggregate to binary data 
# 'cbind' binds the 2 columns 'Coloured' and 'Non-coloured' 
model.form = cbind(Coloured,Total-Coloured)~Housing+Feeding
# Transforming into binary data
lamb.bin = splitbin(model.form,data=lamb)$tab[,-1]
head(lamb.bin)

## Two-way analysis of variance logistic model

lamb.glm = glm(Coloured~Housing*Feeding,data=lamb.bin,family=binomial)
coef(lamb.glm)
odds.ratio(lamb.glm)

## One-way analysis of variance logistic model

lamb.glm = glm(Coloured~Feeding,data=lamb.bin,family=binomial)
coef(lamb.glm)
odds.ratio(lamb.glm)
