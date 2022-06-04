
# PCA

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(mvtnorm)
library(ellipse)

# ------------------------------------------------------------------------------------
# 1. Basic EDA
# ------------------------------------------------------------------------------------
# DATA FORMAT
#
#     |----|----|----|----|
#     | X1 | X2 | .. | Xp |
# |---|----|----|----|----|
# | 1 | .. | .. | .. | .. |
# | 2 | .. | .. | .. | .. |
# |   | .. | .. | .. | .. |
# | n | .. | .. | .. | .. |
# |---|----|----|----|----|
#
# Attention!
# Variables have to be continuous!
# If there are categorical data, get rid of them.

# Read the data (and consider only continuous variables!)
data = read.table('tourists.txt', header=T)
data = data[3:10]
head(data)

# Boxplot
boxplot(data, col='gold')

# Scaled boxplot
boxplot(scale(x=data, center=T, scale=F), col='gold')

# Pairs
pairs(data)

# Matplot
matplot(t(data), type='l', labels=FALSE)
boxplot(data, add=T, boxwex=0.1, col = 'red')

# ------------------------------------------------------------------------------------
# 2. Some estimators:
#   * mean
#   * covariance
#   * correlation
#   * generalized variance
#   * total variance
# ------------------------------------------------------------------------------------
# Mean
M = sapply(data, mean)
M

# Covariance
S = cov(data)
image(S)
round(S, digits=2)

# Correlation
R = cor(data)
round(R, digits=2)

# Generalized variance
var.gen = det(S)
var.gen

# Total variance
var.tot =sum(diag(S))
var.tot

# ------------------------------------------------------------------------------------
# 3. PCA on original data
# ------------------------------------------------------------------------------------
pc.data = princomp(data, scores=T)
pc.data
summary(pc.data)
# Looking at the cumulative proportion:
# is it possible to perform dimensionality reduction?

# Proportion of variance explained by each PC
pc.data$sd^2/sum(pc.data$sd^2)

# Cumulative proportion
cumsum(pc.data$sd^2)/sum(pc.data$sd^2)

# Explained variance plot
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.data, las=2, main='Principal components')
barplot(sapply(data,sd)^2, las=2, main='Original Variables', ylab='Variances')
plot(cumsum(pc.data$sd^2)/sum(pc.data$sd^2), type='b', axes=F,
     xlab='number of components', ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(data),labels=1:ncol(data),las=2)

# Loadings
load.data = pc.data$loadings

# Graphical representation of the loadings
par(mfrow=c(2,4))
for( i in 1:dim(data)[2]){
  barplot(load.data[,i], ylim=c(-1,1), main=paste("PC",i))
}

# ------------------------------------------------------------------------------------
# 4. PCA on standardized data
# ------------------------------------------------------------------------------------
data.sd = scale(data)
data.sd = data.frame(data.sd)

pc.data.sd = princomp(data.sd, scores=T)
pc.data.sd
summary(pc.data.sd)

# Proportion of variance explained by each PC
pc.data.sd$sd^2/sum(pc.data.sd$sd^2)

# Cumulative proportion
cumsum(pc.data.sd$sd^2)/sum(pc.data.sd$sd^2)

# Explained variance
layout(matrix(c(2,3,1,3),2,byrow=T))
plot(pc.data.sd, las=2, main='Principal Components', ylim=c(0,7))
abline(h=1, col='blue')
barplot(sapply(data.sd,sd)^2, las=2, main='Original Variables',
        ylim=c(0,7), ylab='Variances')
plot(cumsum(pc.data.sd$sde^2)/sum(pc.data.sd$sde^2), type='b', axes=F,
     xlab='Number of components', ylab='Contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(data.sd),labels=1:ncol(data.sd),las=2)

# Loadings
load.data.sd = pc.data.sd$loadings

# Graphical representation of the loadings
par(mfrow=c(2,4))
for( i in 1:dim(data)[2]){
  barplot(load.data.sd[,i], ylim=c(-1,1), main=paste("PC",i))
}

# ------------------------------------------------------------------------------------
# 5. Scores
#
# If we want to perform it on the standardized data:
#   data = data.sd
#   pc.data = pc.data.sd
# ------------------------------------------------------------------------------------
scores.data = pc.data$scores
head(scores.data)
summary(scores.data)

# Boxplot original data vs. boxplot PC
layout(matrix(c(1,2),2))
boxplot(data, las=2, col='gold', main='Original data')
scores.data = data.frame(scores.data)
boxplot(scores.data, las=2, col='grey',main='Principal components')

# Plot
plot(scores.data[,1:2])
abline(h=0, v=0, lty=2, col='grey')
text(scores.data[,1], scores.data[,2], dimnames(data)[[1]], cex=0.7)

# Plot
biplot(pc.data)

# ------------------------------------------------------------------------------------
# 6. Projections on the space of Principal COmponents
# 
# DATA FORMAT: 
#
#            |-----------|-----------|----|------------|
#            |   day_ 1  |   day_2   | .. |   day_n    |
#            | /feature1 | /feature2 | .. | /feature_n |
# |----------|-----------|-----------|----|------------|
# | Person_1 |     ..    |     ..    | .. |     ..     |
# | Person_2 |     ..    |     ..    | .. |     ..     |
# |----------|-----------|-----------|----|------------|
#
# If it is transposed, then: data = t(data)
# If we want to work with standardized data:
#   data = data.sd
#   load.data = load.data.sd
# ------------------------------------------------------------------------------------
# Numero delle componenti principali:
numero_pc = dim(data)[2]

# 1. "Projection on the space generatebu the k component"
matplot(t(data), type='l', main = 'Data', ylim=range(data))
meanF = colMeans(data)
matplot(meanF, type='l', main = '0 PC', lwd=2, ylim=range(data))
for(i in 1:numero_pc){
  projection <- matrix(meanF, dim(data)[[1]], dim(data)[[2]], byrow=T) +
    scores.data[,i] %*% t(load.data[,i])
  matplot(t(projection), type='l', main = paste(i, 'PC'), ylim=range(data))
  matplot(meanF, type='l', lwd=2, add=T)
}

# 2. "Projection on the space generated by the first k components"
matplot(t(data), type='l', main = 'Data', ylim=range(data))
meanF <- colMeans(data)
matplot(meanF, type='l', main = 'First 0 PCs', lwd=2, ylim=range(data))
projection <- matrix(meanF, dim(data)[[1]], dim(data)[[2]], byrow=T)
for(i in 1:numero_pc){
  projection <- projection + scores.data[,i] %*% t(load.data[,i])
  matplot(t(projection), type='l', main = paste('First', i, 'PCs'), ylim=range(data))
  matplot(meanF, type='l', lwd=2, add=T)
}

# ------------------------------------------------------------------------------------
# 7. Scores for a new element
# ------------------------------------------------------------------------------------
# New element
new_elem = c(13,10,11,13)

# Scores
scores.new_elem = t(pc.data$loadings) %*% (new_elem-colMeans(data))
scores.new_elem

# Plot
plot(scores.data[,1], scores.data[,2], col='gray', xlab='Comp1', ylab='Comp2', pch=19)
points(scores.new_elem[1], scores.new_elem[2], col='red', pch=19)






