
# Test and CR for the mean of a multivariate Gaussian (one population)

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(car)
library(mvtnorm)
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# Testing
#   1.  Formulate the test (and test gaussianity if needed)
#   2.  Compute the test statistic
#   3a. Set alpha and verify if the test statistic belongs to the rejection region
#   3b. Compute the p-value of the test
# ------------------------------------------------------------------------------------
# Test and CR_(1-alpha) for the mean of a multivariate Gaussian (one population)
#     H0: mu == mu0
#     H1: mu != mu0
# ------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------
# 1. Generic p (# features)
# ------------------------------------------------------------------------------------
# Data
data = read.table('stiff.dat', header=T)

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Gaussianity
mcshapiro.test(data)

# If gaussianity is violated we can remove some outliers:
#     d2 = matrix(mahalanobis(data, data.mean, data.cov))
#     data = data[which(d2<7.5),] ## controlla che sia 7.5 ?
#     n = dim(data)[1]
#     p = dim(data)[2]
#     data.mean = sapply(data, mean)
#     data.cov = cov(data)
#     data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0 = c(1900, 1700, 1500, 1720)

# T2 statistics
data.T2 = n*(data.mean - mu0) %*% data.invcov %*% (data.mean - mu0)

# Fisher's quantile
cfr.fisher = ((n-1)*p/(n-p))*qf(1-alpha, p, n-p)

# Test
data.T2 < cfr.fisher # TRUE -> accept H0

# p-value
P = 1 - pf(data.T2*(n-p)/((n-1)*p), p, n-p)
P

# Plot
xx = seq(0,40,by=0.05)
plot(xx, df(xx*(n-p)/((n-1)*p),p,n-p), type="l", lwd=2, main='Density F(p,n-p)',
     xlab='x*(n-p)/((n-1)*p)', ylab='Density')
abline(h=0, v=data.T2*(n-p)/((n-1)*p), col=c('grey','red'), lwd=2, lty=c(2,1))

# ------------------------------------------------------------------------------------
# 2. Number of features (p) = 2
# ------------------------------------------------------------------------------------
# Data
data = read.table('stiff.dat', header=T)
data = data[1:2]

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0  = c(2000,2000)

# Rejection region centered in mu0 (outside the blue ellipse)
plot(data, asp=1)
ellipse(mu0, shape=data.cov/n, sqrt(cfr.fisher), col='blue', lty=2, center.pch=19)

# Sample mean on the plot: red point
points(data.mean[1], data.mean[2], pch=16, col='red', cex=1.5)

# Confidence region (centered in data.mean) level 100*(1-alpha)%
# { x in R^2 s.t. n*(data.mean-x)' %*% (data.cov)^-1 %*% (data.mean-x) < cfr.fisher }
ellipse(data.mean, data.cov/n, sqrt(cfr.fisher), col='red', lty=2, lwd=2, center.cex=1)

# Confidence region with radius as the quantile of order 1-pval
ellipse(data.mean, data.cov/n, sqrt((n-1)*p/(n-p)*qf(1-as.numeric(P),p,n-p)),
        lty=1,col='dark grey',lwd=2)

# Comment
# * If the sample mean is inside the blue ellipse we accept H0.
#   Meaning: the sample mean is not in the rejection region.
# * If mu0 (mean under H0) is in the confidence region of level 1-alpha
#   then we do not reject H0 at level alpha.
# * The confidece region of level 1-alpha contains all the mu0 that we would
#   accept at level alpha.
# * By def. the confidence region of level 1-alpha produces ellipsoidal regions that
#   contain the true mean 100(1-alpha)% of the times. If H0 is true (i.e mu0 is the
#   true mean), those ellipsoidal regions will contain mu0 100(1-alpha)% of the times.

# ------------------------------------------------------------------------------------
# 3. Intervals among directions (p=2)
# ------------------------------------------------------------------------------------
# Data
data = read.table('stiff.dat', header=T)
data = data[1:2]

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0  = c(2000,2000)

# --------------------------------------
# 1. Among the axes of the coordinates
# --------------------------------------
T2 = cbind(inf    = data.mean - sqrt(cfr.fisher*diag(data.cov)/n),
           center = data.mean,
           sup    = data.mean + sqrt(cfr.fisher*diag(data.cov)/n))
T2

# Plot (the red square is the two CI)
plot(data, asp = 1, main='Confidence and rejection regions')
ellipse(mu0, shape=data.cov/n, sqrt(cfr.fisher), col='blue', lty=2, center.pch=16)
points(data.mean[1], data.mean[2], pch = 16, col='red', cex=1.5)
ellipse(data.mean, shape=data.cov/n, sqrt(cfr.fisher), col='red', lty=2, center.pch=16)
rect(T2[1,1],T2[2,1],T2[1,3],T2[2,3], border='red', lwd=2)

# --------------------------------------
# 2. Among the worst direction
#    Direction along which the T2
#    statistics (univariate)
#    is maximized
# --------------------------------------
data.T2 # maximum T2

worst = solve(cov(data)) %*% (data.mean - mu0)
worst = worst/sqrt(sum(worst^2))
worst

theta.worst = atan(worst[2]/worst[1])+pi
theta.worst

IC.worst <- c(data.mean %*% worst - sqrt(cfr.fisher*(t(worst)%*%cov(data)%*%worst)/n),
              data.mean %*% worst,
              data.mean %*% worst + sqrt(cfr.fisher*(t(worst)%*%cov(data)%*%worst)/n) )
IC.worst

# Projecting mu0 on the worst direction
mu0 %*% worst
(IC.worst[1] < mu0%*%worst) & (mu0%*%worst < IC.worst[2])

# Extremes of IC.worst in the coordinate system (x,y):
x.min = IC.worst[1]*worst
x.max = IC.worst[3]*worst
m1.ort = -worst[1]/worst[2]
q.min.ort = x.min[2] - m1.ort*x.min[1]
q.max.ort = x.max[2] - m1.ort*x.max[1]
abline(q.min.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
abline(q.max.ort, m1.ort, col='forestgreen', lty=2,lwd=1)
m1 = worst[2]/worst[1] # worst direction
abline(0, m1, col='grey35')
segments(x.min[1],x.min[2],x.max[1],x.max[2],lty=1,lwd=2, col='forestgreen')

# Comment
# If we reject the global test (H0), even if the single CI would contain (each)
# the two means, nothing would change, it is not a contraddiction. It means that we
# are rejecting H0 at least in one direction, not necessarily the directions of the coords.

# Bonferroni with a generic k
k = p
cfr.t = qt(1-alpha/(2*k), n-1)
Bf = cbind(inf    = data.mean - cfr.t*sqrt(diag(data.cov)/n),
           center = data.mean,
           sup    = data.mean + cfr.t*sqrt(diag(data.cov)/n))
Bf

rect(Bf[1,1],Bf[2,1],Bf[1,3],Bf[2,3], border='orange', lwd=2)
legend('topleft', c('Rej. Reg.', 'Conf. Reg','T2-sim', 'Bonferroni'),
       col=c('blue','red','red','orange'),lty=c(2,2,1,1),lwd=2)

# Comment
# Case in which we reject H0 global but both components of the mean are inside 
# the single confidence intervals: it may be that, even with Bonferroni's correction
# the single CI accept.

# ------------------------------------------------------------------------------------
# 4. Intervals among directions (p>2)
# ------------------------------------------------------------------------------------
# Confidence regions for the mean of level 100(1-alpha)%
# * we want the CR for the mean (ellipsoidal region)
#   { x in R^4 t.c. n*(data.mean-x)' %*% data.invcov %*% (data.mean-m) < cfr.fisher }
# * characterization of the region: centre, direction of the principal axes, length of the axes

# Data
data = read.table('stiff.dat', header=T)

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0 = c(1900, 1700, 1500, 1720)

# Centre
data.mean

# Direction of the semi-axes of the ellipse
r = sqrt(cfr.fisher)
r*sqrt(eigen(data.cov/n)$values)

# We plot the projections of the ellipsoid in some directions of interest
# (e.g the x and y coordinates)

# We plot simultaneous T2 confidence intervals in each direction of interest
# (with global coverage alpha)
T2 = cbind(inf    = data.mean - sqrt(cfr.fisher*diag(data.cov)/n),
           center = data.mean,
           sup    = data.mean + sqrt(cfr.fisher*diag(data.cov)/n))
T2

matplot(1:p, 1:p, pch='', ylim=range(T2), xlab='Variables', ylab='T2 for a component',
        main='Simultaneous T2 conf. int. for the components')
for(i in 1:p) segments(i,T2[i,1],i,T2[i,3],lwd=3,col=i)
points(1:p, T2[,2], pch=16, col=1:4)

# Is mu0 inside the rectangular region?
points(1:p, mu0, lwd=3, col='orange', pch=19)

# Bonferroni
# (with global level 100(1-alpha)%)
k = p
cfr.t = qt(1 - alpha/(k*2), n-1)
Bf = cbind(inf    = data.mean - cfr.t*sqrt(diag(data.cov)/n),
           center = data.mean,
           sup    = data.mean + cfr.t*sqrt(diag(data.cov)/n))
Bf

matplot(1:p, 1:p, pch='', ylim=range(T2), xlab='Variables',
        ylab='Confidence intervals along a component',main='Confidence intervals')
for(i in 1:p) segments(i,T2[i,1],i,T2[i,3],lwd=2,col='grey35', lty=3)
points(1:p, T2[,1], pch='-', col='grey35')
points(1:p, T2[,3], pch='-', col='grey35')
for(i in 1:p) segments(i,Bf[i,1],i,Bf[i,3],lwd=2,col=i)
points(1:p, Bf[,2], pch=16, col=1:p)
points(1:p, Bf[,1], pch='-', col=1:p)
points(1:p, Bf[,3], pch='-', col=1:p)

# Is mu0 inside the Bonferroni confidence region?
# We add it to the plot
points(1:p, mu0, lwd=3, col='orange', pch=19)

# ------------------------------------------------------------------------------------
# 5. Asymptotic test on the mean (generic p)
# (no gaussianity of the data but n must be high)
# ------------------------------------------------------------------------------------
# Data
data = read.table('stiff.dat', header=T)

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0 = c(1900, 1700, 1500, 1720)

# T2 statistics
data.T2A = n*(data.mean - mu0) %*% data.invcov %*% (data.mean - mu0)

# Radius of the ellipsoid
cfr.chisq = qchisq(1-alpha, p) # Chi-square's quantile

# Test
data.T2A < cfr.chisq # TRUE -> accept H0

# p-value
PA = 1 - pchisq(data.T2A, p)
PA

# Plot
curve(dchisq(x, df = p), from = 0, to = 40, main = 'Chi-Square Distribution',
      ylab = 'Density', lwd = 2)
abline(h=0, v=data.T2A, col=c('grey','red'), lwd=2, lty=c(2,1))

# ------------------------------------------------------------------------------------
# 6. Comparison of the rejection regions (p=2)
#    (Fisher's vs. Chi-square's)
# ------------------------------------------------------------------------------------
# Data
data = read.table('stiff.dat', header=T)
data = data[1:2]

# Dimensions
n = dim(data)[1]
p = dim(data)[2]

# Estimators
data.mean   = sapply(data, mean)
data.cov    = cov(data)
data.invcov = solve(data.cov)

# Choose alpha
alpha = 0.01

# Mean to test
mu0  = c(2000,2000)

# Plot
plot(data, asp = 1,main='Comparison rejection regions')
ellipse(mu0, shape=data.cov/n, sqrt(cfr.fisher), col = 'blue',
        lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
ellipse(mu0, data.cov/n, sqrt(cfr.chisq), col = 'lightblue',
        lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
points(mu0[1], mu0[2], pch = 4, cex = 1.5, lwd = 2, col ='lightblue')
legend('topleft', c('Exact', 'Asymptotic'),col=c('blue','lightblue'),lty=c(1),lwd=2)

# Comparison of the confidence regions
plot(data, asp = 1,main='Comparison of confidence regions')
ellipse(data.mean, data.cov/n, sqrt(cfr.fisher), col = 'red',
        lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
ellipse(data.mean, data.cov/n, sqrt(cfr.chisq), col = 'orange',
        lty = 1, center.pch = 4, center.cex=1.5, lwd=2)
points(data.mean[1], data.mean[2], pch = 4, cex = 1.5, lwd = 2, col ='orange')
legend('topleft', c('Exact', 'Asymptotic'),col=c('red','orange'),lty=c(1),lwd=2)
