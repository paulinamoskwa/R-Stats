
# Functional Data Analysis (FDA)

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(fda)
library(fields)
library(fdakma)

# ------------------------------------------------------------------------------------
# 1. Smoothing
# ------------------------------------------------------------------------------------
# DATA FORMAT
#
#          |------------|------------|----|------------|
#          | location_1 | location_2 | .. | location_n |
# |--------|------------|------------|----|------------|
# | time_1 |     ..     |     ..     | .. |     ..     |
# | time_2 |     ..     |     ..     | .. |     ..     | 
# |   ..   |     ..     |     ..     | .. |     ..     |
# |--------|------------|------------|----|------------|

# If it's not like that (but transposed): data = t(data)

data = CanadianWeather$dailyAv[,,1]
head(data)
dim(data)
n = dim(data)[2]
matplot(data, type='l', xlab='time', ylab='feat')

# ----------------------------------------------
# Fourier basis (periodic)
# ----------------------------------------------
time = 1:dim(data)[1]
nbasis = 10
basis.1 = create.fourier.basis(rangeval = c(0,dim(data)[1]), nbasis = nbasis)
plot(basis.1)
data.fd.1 = Data2fd(y = data, argvals = time, basisobj = basis.1)
plot.fd(data.fd.1)

# "Report the first 3 coefficients of St. Johns"
as.numeric(data.fd.1$coefs[1:3, 'St. Johns'])

# ----------------------------------------------
# Least squares basis (without penalization)
# ----------------------------------------------
time = 1:dim(data)[1]
m = 5 # spline order
degree = m-1 # spline degree
nbasis = 9 # number of basis
basis.2 = create.bspline.basis(rangeval = c(0,dim(data)[1]), nbasis = nbasis)
plot(basis.2)
data.fd.2 = Data2fd(y = data, argvals = time, basisobj = basis.2)
plot.fd(data.fd.2)

# ----------------------------------------------
# Least squares basis (with penalization)
# ----------------------------------------------
time = 1:dim(data)[1]
m = 5 # spline order
degree = m-1 # spline degree
NT = length(time)
breaks = time[((0:floor(NT/2))*2)+1]
basis.3 = create.bspline.basis(breaks, norder = m)
plot(basis.3)
data.fd.3 = Data2fd(y = data, argvals = time, basisobj = basis.3)
plot.fd(data.fd.3)

# ----------------------------------------------
# For one curve
# ----------------------------------------------
# curve: data[,1]
# basis: basis.1

# Smooth curve
basismat0 = eval.basis(time, basis.1)
Xsp0 = basismat0 %*% lsfit(basismat0, data[,1], intercept=F)$coef

plot(time, data[,1])
points(time, Xsp0, type='l', col='blue', lwd=2)

# First derivative
#   finite differences
NT = dim(data)[2]
rappincX1 = (data[3:NT,1]-data[1:(NT-2),1])/(time[3:NT]-time[1:(NT-2)])
basismat1 = eval.basis(time, basis.1, Lfdobj=1)
Xsp1 = basismat1 %*% lsfit(basismat0, data[,1], intercept=F)$coef
plot(time[2:(NT-1)], rappincX1, xlab='t', ylab='first derivative', type='l')
points(time, Xsp1, type='l', col='orange', lwd=3)

# Second derivative
#   finite differences
rappincX2 = ((data[3:NT,1]-data[2:(NT-1),1])/(time[3:NT]-time[2:(NT-1)])-
               (data[2:(NT-1)]-data[1:(NT-2)])/(time[2:(NT-1)]-time[1:(NT-2)]))*
  2/(time[3:(NT)]-time[1:(NT-2)])
basismat2 = eval.basis(time, basis.1, Lfdobj=2)
Xsp2 = basismat2 %*% lsfit(basismat0, data[,1], intercept=F)$coef
plot(time[2:(NT-1)], rappincX2, xlab="t", ylab="second derivative", type="l")
points(time, Xsp2, type='l', col="orange", lwd=3)

# ----------------------------------------------
# Approximate pointwise confidence intervals (one at the time!)
# ----------------------------------------------
# we can estimate the variance of x(t) as: sigma^2*diag[phi*(phi'phi)^{-1}(phi)']
S = basismat0 %*% solve(t(basismat0) %*% basismat0) %*% t(basismat0) #projector
sigmahat = sqrt(sum((data[,1]-data[,1])^2)/(NT-nbasis)) #estimate of sigma
lb = Xsp0 - qnorm(0.975) * sigmahat * sqrt(diag(S))
ub = Xsp0 + qnorm(0.975) * sigmahat * sqrt(diag(S))
plot( time, Xsp0, type="l", col="blue", lwd=2, ylab="")
points(time, lb, type="l", col="red", lwd=2, lty="dashed",)
points(time, ub, type="l", col="red", lwd=2, lty="dashed")

# ----------------------------------------------
# Choosing the number of basis by generalized cross-validation
# ----------------------------------------------
# basis: bspline (m=5)
nbasis = 6:150
gcv = numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis = create.bspline.basis(c(0,dim(data)[1]), nbasis[i], m)
  gcv[i] = smooth.basis(time, data[,1], basis)$gcv
}
plot(nbasis, gcv)
nbasis[which.min(gcv)]

# ------------------------------------------------------------------------------------
# 2. Mean and Covariance
# ------------------------------------------------------------------------------------
# Mean
plot.fd(data.fd.1)
lines(mean.fd(data.fd.1), lwd=3)
plot.fd(data.fd.2)
lines(mean.fd(data.fd.2), lwd=3)
plot.fd(data.fd.3)
lines(mean.fd(data.fd.3), lwd=3)

# Covariance
eval.1 = eval.fd(time, data.fd.1)
image.plot(cov(t(eval.1)))
eval.2 = eval.fd(time, data.fd.2)
image.plot(cov(t(eval.2)))
eval.3 = eval.fd(time, data.fd.3)
image.plot(cov(t(eval.3)))

# ------------------------------------------------------------------------------------
# 3. PCA
# ------------------------------------------------------------------------------------
pca.data = pca.fd(data.fd.1, nharm=5, centerfns=T) # nharm = number of PC's

# PCA compute all the pc's, but only n-1 are not null
plot(pca.data$values, xlab='j', ylab='Eigenvalues')
plot(pca.data$values[1:n], xlab='j', ylab='Eigenvalues')
plot(cumsum(pca.data$values)[1:n]/sum(pca.data$values), xlab='j', ylab='CPV')

# Explained variance
pca.data$varprop

# First PC
plot(pca.data$harmonics[1,], col=1, ylab='FPC1', ylim=c(-0.1,0.08))

# Second PC
plot(pca.data$harmonics[2,], col=2, ylab='FPC2', ylim=c(-0.1,0.08))

# Plot of FPCs as perturbation of the mean
par(mfrow=c(1,2))
plot.pca.fd(pca.data)

# Scatterplot of the scores
par(mfrow=c(1,2))
plot(pca.data$scores[,1], pca.data$scores[,2], xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
points(pca.data$scores[n,1], pca.data$scores[n,2],col=2, lwd=4)
plot(pca.data$scores[,1], pca.data$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca.data$scores[,1], pca.data$scores[,2], dimnames(data)[[2]], cex=1)

# Outliers
head(data)
par(mfrow=c(1,1))
matplot(eval.1, type='l')
lines(eval.1[,35], lwd=4, col=2)

# ----------------------------------------------
# Scores with 3 FPCs
# (commented points if the 12-th is an outlier)
# ----------------------------------------------
layout(cbind(1,2,3))
pca_L = pca.data
plot(pca_L$scores[,1],pca_L$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
# points(pca_L$scores[12,1],pca_L$scores[12,2],col=2, lwd=4)
plot(pca_L$scores[,1],pca_L$scores[,3],xlab="Scores FPC1",ylab="Scores FPC3",lwd=2)
# points(pca_L$scores[12,1],pca_L$scores[12,3],col=2, lwd=4)
plot(pca_L$scores[,2],pca_L$scores[,3],xlab="Scores FPC2",ylab="Scores FPC3",lwd=2)
# points(pca_L$scores[12,2],pca_L$scores[12,3],col=2, lwd=4)