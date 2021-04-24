
setwd("C:/Users/David/Dropbox/DS20/Session4")

# Paired comparisons

## Toy data

toydta = data.frame(Observer=rep(c("J1","J2","J3"),c(2,2,2)),
                    Product=rep(c("A","B"),3),Rating=c(2,4,3,6,6,8))
toydta

## One-way analysis of variance F-test

toydta.lm1 = lm(Rating~Product,data=toydta)
anova(toydta.lm1)

## Two-way analysis of variance F-test

toydta.lm2 = lm(Rating~Product+Observer,data=toydta)

### Extracts the least-squares estimation of coefficients
coef(toydta.lm2)

### Extracts the estimated residual standard deviation
summary(toydta.lm2)$sigma

### Reminds the estimated residual standard deviation of the 1-way ANOVA model
summary(toydta.lm1)$sigma

### Analysis of variance table of the 2-way model
anova(toydta.lm2)

### Reminds the analysis of variance table of the 1-way model
anova(toydta.lm1)

## Equivalent paired t-test
t.test(Rating~Product,data=toydta,paired=TRUE)

# Import data

dta = read.table("pig.txt",stringsAsFactors=TRUE)
str(dta)  # Provides a columnwise overview of the data table

# Scatterplot of LMP against back fat depth

plot(LMP~BFAT,data=dta,bty="l",xlab="Backfat depth (mm)",
              ylab="LMP",cex.lab=1.25,pch=16,
              main="Effect of backfat depth on Lean Meat Percentage")
grid()

# Correlation coefficient

## Scaled variables
bfat.scaled = scale(dta$BFAT)[,1]
lmp.scaled = scale(dta$LMP)[,1]

## Check that the scaled series are centered
mean(bfat.scaled);mean(lmp.scaled)

## Check that the standard deviation of the scaled series are 1
sd(bfat.scaled);sd(lmp.scaled)

## Identification of 'extreme' backfat values 
extreme.bfat = which(abs(bfat.scaled)>1.96)
# Display of the 'outliers'
data.frame(WHICH=extreme.bfat,BFAT=bfat.scaled[extreme.bfat],
           LMP=lmp.scaled[extreme.bfat])

# Scatterplot with scaled values
plot(bfat.scaled,lmp.scaled,bty="l",xlab="Scaled backfat depth (mm)",
     ylab="Scaled LMP",cex.lab=1.25,pch=16,
     main="Effect of backfat depth on Lean Meat Percentage")
points(bfat.scaled[extreme.bfat],lmp.scaled[extreme.bfat],col="orange",pch=16)
grid()

## Correlation coefficient between backfat depth and LMP
cor(dta$BFAT,dta$LMP)

# Least-squares fit

## Least-squares fit of the regression model
lmp.lm = lm(LMP~BFAT,data=dta)

## Extract estimated coefficients 
beta = coef(lmp.lm)
beta

## Adds the regression line on the scatterplot
plot(LMP~BFAT,data=dta,bty="l",xlab="Backfat depth (mm)",
              ylab="LMP",cex.lab=1.25,pch=16,
              main="Effect of backfat depth on Lean Meat Percentage")
abline(beta,lwd=2,col="blue")
grid()

## Adds a legend to the plot
legend("bottomleft",lwd=2,bty="n",col=c("blue"),
       legend=c("Least-squares linear fit"))

## Residual standard deviation
summary(lmp.lm)$sigma

