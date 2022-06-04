
# One-way ANOVA (p=1, g>2)

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------------
data = chickwts

# ------------------------------------------------------------------------------------
# 2. Exploration
# ------------------------------------------------------------------------------------
dim(data)
summary(data)
colnames(data) = c('measure', 'label')
head(data)
g = length(levels(data$label))
attach(data)
plot(label, measure, xlab='treat', ylab='measures')

# ------------------------------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------------------------------
# measure.ij = mu + tau.i + eps.ij eps.ij ~ N(0, sigma^2)
#
# H0: tau.1 = tau.2 = .. = tau.g = 0
# H1: there is a tau.j != 0

# Boxplots
par(mfrow=c(1,2))
barplot(rep(mean(measure),g), names.arg=levels(label), ylim=c(0,max(measure)),
        las=2, col='grey85', main='Model under H0')
barplot(tapply(measure, label, mean), names.arg=levels(label), ylim=c(0,max(measure)),
        las=2, col=rainbow(6),main='Model under H1')

# Check the assumptions
n     = length(label) # total number of obs.
ng    = table(label)  # number of obs. in each group
treat = levels(label) # levels of the treatment
g     = length(treat) # number of levels (i.e., of groups)

# 1. Gaussianity of the groups
Ps = 0*(1:g)
for(i in 1:g){
  Ps[i] = shapiro.test(measure[ label==treat[i]])$p
}
Ps

# 2. Same covariance structure (= same sigma^2)
Var = 0*(1:g)
for(i in 1:g){
  Var[i] = var(measure[label==treat[i]])
}
Var

# test of homogeneity of variances
# H0: sigma.1 = sigma.2 = .. = sigma.g
# H1: there exist i,j s.t. sigma.i!=sigma.j
bartlett.test(measure, label)

# ------------------------------------------------------------------------------------
# 4. One-way ANOVA
# ------------------------------------------------------------------------------------
fit = aov(measure ~ label)
summary(fit)

# How to read the summary:
#              Df    Sum Sq    Mean Sq         F value      Pr(>F)
#   treat     (g-1)  SStreat   SStreat/(g-1)   Fstatistic   p-value [H0: tau.i=0 for every i]
#   Residuals (n-g)  SSres     SSres/(n-g)
#   --
#   Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# ------------------------------------------------------------------------------------
# 5. Which supplement is responsible for this?
# ------------------------------------------------------------------------------------
# To see this, we need to do g*(g-1)/2 comparisons.

## ------------------
## Method 1
## Bonferroni
## ------------------
k      = g*(g-1)/2
alpha  = 0.05
Mediag = tapply(measure, label, mean)
SSres  = sum(residuals(fit)^2)
S      = SSres/(n-g)

# CI for all the differences
ICrange = NULL
for(i in 1:(g-1)){
  for(j in (i+1):g){
    print(paste(treat[i],"-",treat[j]))
    print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) *
                         sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) *
                         sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
    ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) *
                                         sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) *
                                         sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }}

# Comment:
#   If for a CI there is NOT 0, then there is evidence that the treatment had an effect
#   in the two groups.

# Plot
par(mfrow=c(1,2))
plot(label, measure, xlab='treat', ylab='measure', col = rainbow(6), las=2)
h = 1
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat',
     ylab='Conf. Int. tau weight', main='Univariate Conf. Int. - Bonf. corrected')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    ind = (i-1)*g-i*(i-1)/2+(j-i)
    lines(c(h,h), c(ICrange[ind,1],ICrange[ind,2]), col='grey55');
    points(h, Mediag[i]-Mediag[j], pch=16, col='grey55');
    points(h, ICrange[ind,1], col=rainbow(g)[j], pch=16);
    points(h, ICrange[ind,2], col=rainbow(g)[i], pch=16);
    h <- h+1
  }}
abline(h=0)

## ------------------
## Method 2
## k one-at-the-time CI
## (without Bonferroni correction)
## ------------------
# In this case each CI has the alpha power, but not together.
# So, we have to change criterium to control the univariate rejection (multiple testing).

Auni = matrix(0,g,g)
for(i in 1:g) {
  for(j in i:g) {
    Auni[i,j] = Mediag[i]-Mediag[j]+qt(1-alpha/2,n-g)*sqrt(S*(1/ng[i]+1/ng[j]))}
  for(j in 1:i) {
    Auni[i,j] = Mediag[j]-Mediag[i]-qt(1-alpha/2,n-g)*sqrt(S*(1/ng[i]+1/ng[j]))}
  Auni[i,i] = 0
}

par(mfrow=c(1,2))
h = 1
plot(c(1,g*(g-1)/2),range(Auni), pch='', xlab='pairs treat',
     ylab='CI delta measure', main='Univariate Conf. Int. - 1-at-the-time', col='grey55')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {lines (c(h,h), c(Auni[i,j],Auni[j,i]));
    points(h, as.numeric(Mediag[i]-Mediag[j]), pch=16, col='grey55');
    points(h, Auni[i,j], col=rainbow(g)[i], pch=16);
    points(h, Auni[j,i], col=rainbow(g)[j], pch=16);
    h <- h+1
  }}
abline(h=0)

# For the comparison we add: 
h = 1
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat',
     ylab='Conf. Int. tau weight', main='Univariate Conf. Int. - Bonf. corrected')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    ind = (i-1)*g-i*(i-1)/2+(j-i)
    lines(c(h,h), c(ICrange[ind,1],ICrange[ind,2]), col='grey55');
    points(h, Mediag[i]-Mediag[j], pch=16, col='grey55');
    points(h, ICrange[ind,1], col=rainbow(g)[j], pch=16);
    points(h, ICrange[ind,2], col=rainbow(g)[i], pch=16);
    h <- h+1
  }}
abline(h=0)

# We compute the p-values of the univariate tests
# Matrix of tests for the difference between all the pairs
P <- matrix(0,g,g)
for(i in 1:g) {
  for(j in i:g) {
    P[i,j] = (1-pt(abs((Mediag[i]-Mediag[j])/sqrt(S*(1/ng[i]+1/ng[j]))), n-g))*2}
  for(j in 1:i) {
    P[i,j] = (1-pt(abs((Mediag[i]-Mediag[j])/sqrt(S*(1/ng[i]+1/ng[j]))), n-g))*2}
  P[i,i] = 0
}
P

# Vector of p-values
p_values <- c(P[1, 2:6], P[2, 3:6], P[3, 4:6], P[4, 5:6], P[5, 6])
p_values

# Plot of all the p-values (one-at-time, Bonferroni, FDR)
plot(1:(g*(g-1)/2), p_values, ylim=c(0,1), type='b', pch=16, col='grey55',
     xlab='pairs treat', main='P-values')
abline(h=alpha, lty=2)

# Bonferroni correction
p.bonf <- p.adjust(p_values, 'bonf')
lines(1:(g*(g-1)/2), p.bonf, col='blue', pch=16, type='b')

# Correction according to the false discovery rate (Benjamini-Hockberg)
p.fdr <- p.adjust(p_values, 'fdr')
lines(1:(g*(g-1)/2), p.fdr, col='red', pch=16, type='b')
legend('topleft', c('Not corr.', 'Bonf.', 'BH'), col=c('grey55', 'blue', 'red'), pch=16)

# Which did make effect?
which(p.bonf < alpha)
which(p.fdr < alpha)
detach(data)