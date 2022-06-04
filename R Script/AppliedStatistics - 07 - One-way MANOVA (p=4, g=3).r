
# One-way MANOVA (p>1, g>2)
# In this particular example g=3

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------------
data = iris

# ------------------------------------------------------------------------------------
# 2. Exploration
# ------------------------------------------------------------------------------------
head(data)
p = dim(data)[2] - 1
v = 'X1'

# Re-name
for(i in 2:p){
  v = c(v, paste('X', i, sep=''))
}
colnames(data) = c(v, 'label')
head(data)

# Dimensions
n     = length(data$label) # total number of obs.
ng    = table(data$label)  # number of obs. in each group
treat = levels(data$label) # levels of the treatment
g     = length(treat)      # number of levels (i.e., of groups)

data.feats = data[,1:p]

# Exploration
colore = rep(rainbow(p), each=50)
pairs(data.feats, col=colore, pch=16)

# Set some indeces
i1 = which(data$label == treat[1])
i2 = which(data$label == treat[2])
i3 = which(data$label == treat[3])

# Plot: different panels -> different group
par(mfrow=c(1,3))
boxplot(data.feats[i1,], main='group 1', ylim=c(0,8), col = rainbow(p))
boxplot(data.feats[i2,], main='group 2', ylim=c(0,8), col = rainbow(p))
boxplot(data.feats[i3,], main='group 3', ylim=c(0,8), col = rainbow(p))

## Plot: different panels -> different feature
par(mfrow=c(1,4))
boxplot(data.feats[,1]~data$label, main='X1', ylim=c(0,8), col = rainbow(3))
boxplot(data.feats[,2]~data$label, main='X2', ylim=c(0,8), col = rainbow(3))
boxplot(data.feats[,3]~data$label, main='X3', ylim=c(0,8), col = rainbow(3))
boxplot(data.feats[,4]~data$label, main='X4', ylim=c(0,8), col = rainbow(3))

# ------------------------------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------------------------------
# measure.ij = mu + tau.i + eps.ij eps.ij ~ N(0, sigma^2) (in R^p)
#
# H0: tau.1 = tau.2 = .. = tau.g = 0
# H1: there is a tau.j != 0

n1 = length(i1)
n2 = length(i2)
n3 = length(i3)
n  = n1+n2+n3

# Check the assumptions
# 1. normality
Ps = NULL
for(i in 1:g){
  Ps = c(Ps, mcshapiro.test(data[get(paste('i',i, sep='')),1:p])$p)
}
Ps

# 2. same covariance structure (homoschedasticity)
S  = cov(data.feats)
S1 = cov(data.feats[i1,])
S2 = cov(data.feats[i2,])
S3 = cov(data.feats[i3,])

# Qualitatively:
round(S1,digits=1)
round(S2,digits=1)
round(S3,digits=1)
par(mfrow=c(1,3))
image(S1, col=heat.colors(100),main='Cov. S1', asp=1, axes = FALSE,
      breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S2, col=heat.colors(100),main='Cov. S2', asp=1, axes = FALSE,
      breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))
image(S3, col=heat.colors(100),main='Cov. S3', asp=1, axes = FALSE,
      breaks = quantile(rbind(S1,S2,S3), (0:100)/100, na.rm=TRUE))

# Comment
# If they're different, note it down and try to go on.

# ------------------------------------------------------------------------------------
# 4. One-way MANOVA
# ------------------------------------------------------------------------------------
fit = manova(as.matrix(data.feats) ~ data$label)
summary.manova(fit, test="Wilks")
# Pr(>F) = p-value of H0 vs. H1
# If it is very small -> H1 -> the treatment was effective

# Comment
# If p<=2 and g<=3 we have an exact test (with Wilks)

# ------------------------------------------------------------------------------------
# 5. Which supplement is responsible?
# (In the case we accept H1)
# ------------------------------------------------------------------------------------
# First of all:
#   Let's see on which variables the group has an effect.
#   Via ANOVA: for each feature we perform an ANOVA to see if the belonging to
#   the group has an effect on the mean of the variables.
summary.aov(fit)

# Comment
# Pr(>F) = p-value small -> the group has an influence on that X_k
# This analysis does NOT say either which groups differ nor which are the variables
# for which the groups differ. 

## Bonferroni CI
#  We want to know the level (of the labels) that introduces the difference.
#  We have to create g*(g-1)/2 intervals.

alpha = 0.05
k     = p*g*(g-1)/2
qT    = qt(1-alpha/(2*k), n-g)
W     = summary.manova(fit)$SS$Residuals

m  = sapply(data.feats,mean)      # estimates mu
m1 = sapply(data.feats[i1,],mean) # estimates mu.1=mu+tau.1
m2 = sapply(data.feats[i2,],mean) # estimates mu.2=mu+tau.2
m3 = sapply(data.feats[i3,],mean) # estimates mu.3=mu+tau.3

inf12 = m1-m2 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
sup12 = m1-m2 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n2) )
inf13 = m1-m3 - qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
sup13 = m1-m3 + qT * sqrt( diag(W)/(n-g) * (1/n1+1/n3) )
inf23 = m2-m3 - qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
sup23 = m2-m3 + qT * sqrt( diag(W)/(n-g) * (1/n2+1/n3) )
CI <- list(g1_g2 = cbind(inf12, sup12),
           g1_g3 = cbind(inf13, sup13),
           g2_g3 = cbind(inf23, sup23))
CI

# Now we have a complete frame (intervals for all the components of tau_i)

# Comment
#   From these intervals we can see who is responsable for the change:
#   if the zero is NOT present in a comparison 'gi_gj' in an interval 'Xk'
#   then the variable 'Xk' is influenced by groups 'gi' and 'gj'.

# Plot: different panels -> different features
par(mfrow=c(2,4))
boxplot(data.feats[,1]~data$label, main=paste(v[1]), ylim=c(0,8), col = rainbow(g))
boxplot(data.feats[,2]~data$label, main=paste(v[2]), ylim=c(0,8), col = rainbow(g))
boxplot(data.feats[,3]~data$label, main=paste(v[3]), ylim=c(0,8), col = rainbow(g))
boxplot(data.feats[,4]~data$label, main=paste(v[4]), ylim=c(0,8), col = rainbow(g))

mg = rbind(m1,m2,m3)
sp.name = v
for(k in 1:p){
  plot(c(1,g*(g-1)/2),ylim=c(-4,4), xlim=c(1,3), pch='',
       xlab='pairs treat', ylab=paste('CI tau',k), main=paste('CI tau',sp.name[k]))
  lines (c(1,1), c(CI[[1]][k,1],CI[[1]][k,2]));
  points(1, mg[1,k]-mg[2,k], pch=16);
  points(1, CI[[1]][k,1], col=rainbow(g)[2], pch=16);
  points(1, CI[[1]][k,2], col=rainbow(g)[1], pch=16);
  lines (c(2,2), c(CI[[2]][k,1],CI[[2]][k,2]));
  points(2, mg[1,k]-mg[3,k], pch=16);
  points(2, CI[[2]][k,1], col=rainbow(g)[3], pch=16);
  points(2, CI[[2]][k,2], col=rainbow(g)[1], pch=16);
  lines (c(3,3), c(CI[[3]][k,1],CI[[3]][k,2]));
  points(3, mg[2,k]-mg[3,k], pch=16);
  points(3, CI[[3]][k,1], col=rainbow(g)[3], pch=16);
  points(3, CI[[3]][k,2], col=rainbow(g)[2], pch=16);
  abline(h=0)
}

# Comment
#   If no one contains 0 it means that every group has relevance in every feature.