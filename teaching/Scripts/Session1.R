
setwd("C:/Users/David/Dropbox/DS20/Session1")

### Required packages

require(RcmdrMisc)    # For numSummary

# Import data

dta = read.table("pig.txt",stringsAsFactors=TRUE)
str(dta)  # Provides a columnwise overview of the data table

# Effect of genetic type on BFAT

## Group summaries

numSummary(dta$BFAT,groups=dta$GENET)

## Mean and median

x = 1:4   # x is the series 1,2,3,4
mean(x)   # mean value of the series
median(x) # median value

x[4] = 40 # The 4th value is now 40
mean(x)
median(x)

## Percentiles

x = c(89.12,89.70,89.73,89.85,89.88,89.90,89.92,89.96,90.04,90.07,90.09,
      90.14,90.17,90.19,90.24,90.25,90.38,90.40,90.48,90.57,90.57,91.24)

length(x)
quantile(x)

## Standard deviation

x = dta$BFAT
x-mean(x)
sum(x-mean(x))

sd(x)
sqrt(sum((x-mean(x))^2)/59)

## Boxplot

boxplot(BFAT~GENET,data=dta,col="darkgray",cex.lab=1.25,pch=16,
   main="Distribution of the backfat depth (mm) across genetic types")
grid()

