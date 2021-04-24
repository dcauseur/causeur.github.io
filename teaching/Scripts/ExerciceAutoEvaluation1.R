
# Set working directory

setwd("C:/Users/David/Dropbox/DS20/Session1")

# Import data

coffee = read.table("coffee.txt",stringsAsFactors=TRUE)
str(coffee)

# Convert "Localisation" into a factor

coffee$Localisation = factor(coffee$Localisation)
str(coffee)

