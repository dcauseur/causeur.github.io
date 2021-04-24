
### Set working directory

setwd("c:/users/david/dropbox/sad2020/session5/data")

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

install.packages("ROCR")   # Only if needed
require(ROCR)              # For classification error rates

install.packages("groupdata2")   # Only if needed
require(groupdata2)              # For cross-validation

### Import data

# Import data stored in file maturity.txt
dta = read.table("maturity.txt")
# Convert the Maturity variable into a 3-level factor
dta$Maturity = factor(dta$Maturity,levels=c('1','2','3'))
# Overview of data
str(dta)

### Prediction using a logistic model

## 3-group maturity

# Fit a multinomial logistic regression model for Maturity
maturity.mlogit = multinom(Maturity~.,data=dta[,-7])

# Prediction for apricots 1, 43, 211 within the sample  
probas = predict(maturity.mlogit,newdata=dta[c(1,43,211),],type="probs")
probas
classes = predict(maturity.mlogit,newdata=dta[c(1,43,211),],type="class")
classes

## 2-group maturity

dta12 = subset(dta, Maturity != "3")
# Keeps only rows with Maturity != "3"
dta12 = droplevels(dta12)
# Forgets the unused levels
str(dta12)

# Fit a logistic regression model for Maturity
maturity.logit = glm(Maturity~.,data=dta12[,-7],family=binomial)

# Prediction for apricots 1, 43, 211 within the sample  
probas = predict(maturity.logit,newdata=dta[c(1,43,211),],type="response")

# Prediction for all apricots within the sample  
probas = predict(maturity.logit,type="response")
plot(probas~dta12$Maturity,bty="l",col="coral4",xlab="Maturity",
     ylab="Estimated probability",cex.lab=1.25,
     cex.axis=1.25,cex.main="1.25",pch=16,
     main="Probability of being in maturity group 2")

# Bayes clasification rule
pred.class = ifelse(probas>=0.5,"2","1")
confusion = table(dta12$Maturity,pred.class,dnn=list("Observed","Predicted"))
confusion

# Misclassification error rate
mean(pred.class!=dta12$Maturity)

# Prediction performance criteria
mean(pred.class[dta12$Maturity=="2"]=="2") # TPR
mean(pred.class[dta12$Maturity=="1"]=="1") # TNR
mean(dta12$Maturity[pred.class=="2"]=="2") # PPV
mean(dta12$Maturity[pred.class=="1"]=="1") # NPV

# Create a 'prediction' object to be used by function 'performance'
pred = prediction(predictions=probas,labels=dta12$Maturity)

# Derive the prediction performance indicators for a sequence of decision cutoffs
tpr = performance(pred,measure="tpr")
tnr = performance(pred,measure="tnr")
ppv = performance(pred,measure="ppv")
npv = performance(pred,measure="npv")

# Decision cutoffs
cutoffs = tpr@"x.values"[[1]]
cutoffs[1:10]

# Prediction performance indicators for cutoff close to 0.50
choice = which.min(abs(cutoffs-0.5))
cutoffs[choice]
tpr@"y.values"[[1]][choice]
tnr@"y.values"[[1]][choice]
ppv@"y.values"[[1]][choice]
npv@"y.values"[[1]][choice]

# Plots TPR and FPR against the threshold

plot(tpr,lwd=2,col="blue",ylab="TPR and TNR",xlab="Decision threshold")
plot(tnr,lwd=2,col="orange",add=TRUE)
legend(0.4,0.25,bty="n",lwd=2,col=c("blue","orange"),legend=c("TPR","TNR"))

# Finds the minimal threshold for which TPR>=0.95
choice = min(which(tpr@"y.values"[[1]]>=0.95))
cutoffs[choice]
abline(v=cutoffs[choice],lwd=2,col="darkgray",lty=5)
tpr@"y.values"[[1]][choice] # Corresponding TPR
abline(h=tpr@"y.values"[[1]][choice],lwd=2,col="cadetblue",lty=5)
tnr@"y.values"[[1]][choice] # Corresponding TNR
abline(h=tnr@"y.values"[[1]][choice],lwd=2,col="coral",lty=5)

# ROC curve
perf = performance(pred,measure="tpr",x.measure="fpr")
plot(perf,lwd=2,col="blue")

# Ideal point with FPR=0 and TPR=1
points(0,1,pch=16,col="coral",cex=1.5)
text(0.02,1,"Ideal point",cex=1.25,adj=0,col="coral")

# Worst ROC curve
abline(a=0,b=1,lwd=2,col="darkgray")
text(0.6,0.6-0.02,"Worst ROC curve",adj=0,col="darkgray",cex=1.25)

# Area under the ROC curve
performance(pred,measure="auc")@"y.values"[[1]]

# Best compromise TPR-FPR (FPR=1-TNR)
fpr = performance(pred,measure="fpr")
choice = which.min(fpr@"y.values"[[1]]^2+(1-tpr@"y.values"[[1]])^2)
cutoffs[choice]

tpr@"y.values"[[1]][choice] # Corresponding TPR
fpr@"y.values"[[1]][choice] # Corresponding FPR

points(fpr@"y.values"[[1]][choice],tpr@"y.values"[[1]][choice],cex=1.5,pch=16,
       col="cadetblue")
text(fpr@"y.values"[[1]][choice]+0.02,tpr@"y.values"[[1]][choice]-0.02,
     "Closest point to the ideal",cex=1.25,adj=0,col="cadetblue")

### Simple cross-validation setup

# First, random permutation of sampling items
n = nrow(dta12)
permutation = sample(1:n)

# Then, split the data in 2/3 for learning, 1/3 for test
learn_ids = permutation[1:round((2/3)*n)]

learn = dta12[learn_ids,-7] 
test = dta12[-learn_ids,-7]   

# Fit the whole classification rule on the learning sample
maturity.select = bestglm(Xy=learn,family=binomial,
                          method="exhaustive",IC="AIC")
# Extract AIC for best models with k explanatory variables
AIC = maturity.select$Subsets$AIC
# Identifies column numbers for selected explanatory variables
selected = unlist(maturity.select$Subsets[which.min(AIC),2:6])
# Fit the model with selected axplanatory variables
maturity.logit = glm(Maturity~.,family=binomial,
                     data=learn[,c(selected,TRUE)])
# Calculate CV'd probabilities of being in maturity group 2
cvprobas = predict(maturity.logit,newdata=test,type="response")

# Cross-validated AUC
cvpred = prediction(predictions=cvprobas,labels=test$Maturity)
cvAUC = performance(cvpred,measure="auc")@"y.values"[[1]]
cvAUC

### K-fold cross-validation

# Step 1: segmentation of the dataset in 10 segments
segments = fold(dta12,k=10,cat_col="Maturity")$".folds"
# 10-fold balanced partition of the sample
# segments is a vector giving the segment number of each item
table(dta12$Maturity,segments)

# Step 2: cycling over the segments
cvprobas = rep(0,n) # will contain the cross-validated probabilities
for (k in 1:10) {          # cycling over the 10 segments
  # the kth segment is excluded from the learning sample
  learn = dta12[segments!=k,-7] 
  # the test sample is just the kth segment
  test = dta12[segments==k,-7]   
  # Fit the whole classification rule on the learning sample
  maturity.select = bestglm(Xy=learn,family=binomial,
                            method="exhaustive",IC="AIC")
  # Extract AIC for best models with k explanatory variables
  AIC = maturity.select$Subsets$AIC
  # Identifies column numbers for selected explanatory variables
  selected = unlist(maturity.select$Subsets[which.min(AIC),2:6])
  # Fit the model with selected axplanatory variables
  maturity.logit = glm(Maturity~.,family=binomial,
                       data=learn[,c(selected,TRUE)])
  # Calculate CV'd probabilities of being in maturity group 2
  cvprobas[segments==k] = predict(maturity.logit,newdata=test,type="response")
}

# Step 3: Cross-validated performance criteria

cvpred = prediction(predictions=cvprobas,labels=dta12$Maturity)
cvAUC = performance(cvpred,measure="auc")@"y.values"[[1]]
cvAUC

