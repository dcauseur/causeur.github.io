
setwd("C:/Users/David/Dropbox/DS20/Session2")

# Import data

dta = read.table("pig.txt",header=TRUE,stringsAsFactors=TRUE)
str(dta)  # Provides a columnwise overview of the data table

# Effect of genetic type on BFAT

## Least-squares fit of the 1-way ANOVA model

bfat.lm = lm(BFAT~GENET,data=dta)

### Extract estimated coefficients
coef(bfat.lm)

### Summary extracts many useful statistics from the fitted model
### including the residual standard deviation
summary(bfat.lm)$sigma

## ANOVA F-test

### Null versus nonnull model

bfat.lm0 = lm(BFAT~1,data=dta)
coef(bfat.lm0)

RSS = sum(residuals(bfat.lm)^2)
RSS0 = sum(residuals(bfat.lm0)^2)

### Extracts the R2
summary(bfat.lm)$r.squared

### Extracts the F-test statistics
summary(bfat.lm)$fstatistic

### Displays the complete ANOVA table
anova(bfat.lm)

