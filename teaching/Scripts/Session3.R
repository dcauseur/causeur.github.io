
setwd("C:/Users/David/Dropbox/DS20/Session3")

# Import data

dta = read.table("pig.txt",header=TRUE,stringsAsFactors=TRUE)
str(dta)  # Provides a columnwise overview of the data table

# F-test for the effect of sex on BFAT

bfatsex.lm = lm(BFAT~SEX,data=dta)
anova(bfatsex.lm)

# Equivalent t-test 

## var.equal=TRUE in t.test states that 
## within-group standard deviation are assumed to be equal
## mu=0 states that the difference in means under the null is zero
## mu=0 is the default option
t.test(BFAT~SEX,var.equal=TRUE,data=dta,mu=0)$statistic

## The p-value is provided as an output of t.test
t.test(BFAT~SEX,var.equal=TRUE,data=dta,mu=0)$p.value

# One-sided t-test

## alternative="less" is used for the present one-sided test
t.test(BFAT~SEX,var.equal=TRUE,data=dta,alternative="less")$p.value

# Confidence intervals for the mean difference

## Largest value of |t| for which the null is not rejected
qt(0.975,df=58)
qt(0.025,df=58)

t.test(BFAT~SEX,var.equal=TRUE,data=dta)

# Confidence intervals for mean backfat depths of females 
t.test(dta[dta$SEX=="F","BFAT"]) 

# Power of the t-test

sigma = summary(bfatsex.lm)$sigma # Residual standard deviation

power.t.test(delta=1,n=30,sd=sigma,sig.level=0.05)
## delta is the mean difference to detect at population level
## n = 30 is the group size
## sd is the within-group standard deviation
## sig.level is the type-I error level

power.t.test(delta=5,n=30,sd=sigma,sig.level=0.05)$power

power.t.test(delta=1,power=0.90,sd=sigma,sig.level=0.05)$n

# Pairwise t-tests

pairwise.t.test(x=dta$BFAT,g=dta$GENET,p.adjust.method="none")
   # p-values of the pairwise comparisons

pairwise.t.test(x=dta$BFAT,g=dta$GENET,p.adjust.method="bonferroni")
   # adjusted p-values of the pairwise comparisons

