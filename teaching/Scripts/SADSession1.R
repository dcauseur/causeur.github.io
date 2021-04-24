
### Set working directory

setwd("c:/users/david/dropbox/sad2020/session1/data")

### Covid blood group study

# Contingency table for Wuhan's data
covidblood_wuhan = matrix(c(1188,920,336,1250,670,469,178,458),
                          nrow=2,byrow=TRUE,
                          dimnames=list(c("Controls","Cases"),c("A","B","AB","O")))
covidblood_wuhan

# Total number of recruitments
n = sum(covidblood_wuhan)
n

# Marginal totals for control/case
covid.margin = rowSums(covidblood_wuhan)
covid.margin

# Marginal totals for blood types
blood.margin = colSums(covidblood_wuhan)
blood.margin

# Blood type distribution in sample
blood.margin/n

# Expected counts under independence
tab_star = outer(covid.margin,blood.margin)/n
tab_star

# Same with function chisq.test
chisq.test(covidblood_wuhan)$expected
covidblood_wuhan

# Pearson's chi-square test statistic
# Correct = FALSE to get the original test statistic
chisq.test(covidblood_wuhan,correct=FALSE)$statistic

# p-value
chisq.test(covidblood_wuhan,correct=FALSE)

# Residuals
chisq.test(covidblood_wuhan,correct=FALSE)$residuals

# Goup probabilities of infection
covidblood_wuhan["Cases",]/colSums(covidblood_wuhan)

# Odds for each blood type
odds = covidblood_wuhan["Cases",]/covidblood_wuhan["Controls",]
odds

# Odds-ratio of group O versus group A+B+AB
odds_nonO = (sum(covidblood_wuhan["Cases",c("A","B","AB")])/
                sum(covidblood_wuhan["Controls",c("A","B","AB")])) 
odds["O"]/odds_nonO


