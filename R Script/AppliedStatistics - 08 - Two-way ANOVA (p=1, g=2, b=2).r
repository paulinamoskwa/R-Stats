
# Two-way ANOVA (p=1, g=2, b=2)

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# 1. Create the dataset
# ------------------------------------------------------------------------------------
measure = c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
label_1 = factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
label_2 = factor(c('95','95','98','98','95','95','98','98'))

# ------------------------------------------------------------------------------------
# 2. Exploration
# ------------------------------------------------------------------------------------
# Parameters
g = length(levels(label_1))
b = length(levels(label_2))
n = length(measure)/(g*b)

# Means
M = mean(measure)
M_label_1 = tapply(measure, label_1, mean)
M_label_2 = tapply(measure, label_2, mean)

# !Just for visualization!
label_1_2 = factor(c('Esso95','Esso95','Esso98','Esso98','Shell95',
                     'Shell95','Shell98','Shell98'))
M_label_1_2 = tapply(measure, label_1_2, mean)

# Plot
par(mfrow=c(2,3),las=2)
barplot(rep(M,4), names.arg=levels(label_1_2), ylim=c(0,24), main='No factor')
barplot(rep(M_label_1,each=2), names.arg=levels(label_1_2), ylim=c(0,24),
        col=rep(c('blue','red'),each=2), main='label_1')
barplot(rep(M_label_2,times=2), names.arg=levels(label_1_2), ylim=c(0,24),
        col=rep(c('darkgreen','orange'),times=2), main='label_2')
barplot(c(M_label_1[1]+M_label_2[1]-M, M_label_1[1]+M_label_2[2]-M,
          M_label_1[2]+M_label_2[1]-M, M_label_1[2]+M_label_2[2]-M),
        names.arg=levels(label_1_2), ylim=c(0,24),
        col=rep(c('darkgreen','orange'),times=2), density=rep(10,4), angle=135,
        main='Additive model label_1 + label_2')
barplot(c(M_label_1[1]+M_label_2[1]-M, M_label_1[1]+M_label_2[2]-M,
          M_label_1[2]+M_label_2[1]-M, M_label_1[2]+M_label_2[2]-M),
        names.arg=levels(label_1_2), ylim=c(0,24),
        col=rep(c('blue','red'),each=2), density=rep(10,4), add=T)
barplot(M_label_1_2, names.arg=levels(label_1_2), ylim=c(0,24),
        col=rainbow(5)[2:5], main='Model with Interact. label_1 & label_2.')
plot(label_1_2, measure, col=rainbow(5)[2:5], ylim=c(0,24),xlab='')

# ------------------------------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------------------------------
# Note that: if we have to remove rows -> one at the time

###------------------------------
### Complete model
###------------------------------
# measure.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk eps.ijk ~ N(0, sigma^2)
# i=1,2 (label_1)
# j=1,2 (label_2)

fit.aov2.int = aov(measure ~ label_1 + label_2 + label_1:label_2)
summary.aov(fit.aov2.int)

###------------------------------
### Additive model
###------------------------------
# measure.ijk = mu + tau.i + beta.j + eps.ijk eps.ijk ~ N(0, sigma^2)
# i=1,2 (label_1)
# j=1,2 (label_2)

fit.aov2.ad = aov(measure ~ label_1 + label_2)
summary.aov(fit.aov2.ad)

###------------------------------
### Reduced additive model 
### (ANOVA one-way, b=2 (or g=2))
###------------------------------
# measure.jk = mu + beta.j + eps.jk eps.ijk ~ N(0, sigma^2)
# j=1,2 (label_2)

fit.aov1 = aov(measure ~ label_2)
summary.aov(fit.aov1)

# ------------------------------------------------------------------------------------
# 4. IC for the reduced additive model
# ------------------------------------------------------------------------------------
SSres = sum(residuals(fit.aov1)^2)
IC    = c(diff(M_label_2)-qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) * (1/(n*g)+1/(n*g))),
          diff(M_label_2)+qt(0.95, (n*g-1)*b) * sqrt(SSres/((n*g-1)*b) * (1/(n*g)+1/(n*g))))
names(IC) = c('Inf', 'Sup')
IC # IC for mu(label_2[1]) - mu(label_2[2])

# Note
#   There should NOT be the zero.
# ------------------------------------------------------------------------------------
# 5. (Approximate) check of the assumptions
# ------------------------------------------------------------------------------------
# Check one-by-one and just for the considered labels

# 1) Gaussianity
treat = levels(label_2)
Ps = 0*(1:b)
for(i in 1:b){
  Ps[i] = shapiro.test(measure[label_2==treat[i]])$p
}
Ps

# 2) homogeneity of variances
bartlett.test(measure, label_2)
