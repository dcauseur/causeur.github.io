
# Set working directory

setwd("C:/Users/David/Dropbox/DS20/Session1")

# Import data

coffee = read.table("coffee.txt",stringsAsFactors=TRUE)
str(coffee)

# Convert "Localisation" into a factor

coffee$Localisation = factor(coffee$Localisation)
str(coffee)

# Create a new data frame restricted to sites "1" and "6"

coffee2 = subset(coffee,(Localisation=="1")|(Localisation=="6"))

# Refresh the factors by dropping the unused levels

coffee2 = droplevels(coffee2)
str(coffee2)

