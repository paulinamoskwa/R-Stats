
# Paired Gaussian Data: Test for the Mean

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(car)
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# 1. Data settings
# ------------------------------------------------------------------------------------
data = read.table('effluent.dat')

# DATA FORMAT
#   k columns that will be paired
#
#     |----|----|----|----|
#     | V1 | V2 | V3 | V4 |
# |---|----|----|----|----|
# | 1 |    |    |    |    |
# | 2 |    |    |    |    |
# |   |    |    |    |    |
# |---|----|----|----|----|
#
# We're working on a multidimensional mean composed of k/2 components.
# For example, here we'll have: 
#     mu = (V1-V3, V2-V4)
# and we'll want to check if it equal to some mu0

# Pairs
pairs(data, pch=19)

# Rewriting of the dataset
D = data.frame(X1 = data$V1 - data$V3, X2 = data$V2 - data$V4)
D

# Plot
plot(D, asp=1, pch=19, main='dataset of differences')
abline(h=0, v=0, col='grey')
points(0,0, pch=19, col='grey')
# Do we have enough evidence to say that the grey point is the mean?

# ------------------------------------------------------------------------------------
# 2. Work on the mean
# ------------------------------------------------------------------------------------
data = D

# From here it is sufficient to use the previous code.
# ('Test and CR for the mean of a multivariate Gaussian - One population')