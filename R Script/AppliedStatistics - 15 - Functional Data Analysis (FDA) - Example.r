
# Functional Data Analysis (FDA)

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(fda)

# ------------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------------
# File watertemp.txt contains the mean daily water temperature registered at
# 132 monitoring stations in the Adriatic Sea, during the 365 days of 2017.
# The dataset also report the zone of the measurement (Deep, Medium or Surface water).
data = read.table('watertemp.txt', header=T)
data = t(data[,-366])

# ------------------------------------------------------------------------------------
# 2. Request a)
# ------------------------------------------------------------------------------------
# Perform a smoothing of the data through a projection over a Fourier basis with 45 
# basis elements. Report the first 3 Fourier coefficients obtained at the Stations 1 and 2.
dim(data)
basis.1 = create.fourier.basis(rangeval=c(0,365), nbasis=45)
time = 1:365
data.fd.1 = Data2fd(y = data, argvals=time, basisobj = basis.1)
plot.fd(data.fd.1)
as.numeric(data.fd.1$coefs[1:3, 'Station1'])
as.numeric(data.fd.1$coefs[1:3, 'Station2'])

# ------------------------------------------------------------------------------------
# 3. Request b)
# ------------------------------------------------------------------------------------
# Perform a functional principal component analysis of the smoothed data obtained
# at the previous point. Report the variance explained along the first 5 functional
# principal components, a qualitative plot of the first 3 eigenfunctions
# and the screeplot. Interpret the principal components.

# PCA
pca_data = pca.fd(data.fd.1, n=5, centerfns=TRUE)

# Scree plot
plot(pca_data$values)

# Variance explained along the first 5 FPC
pca_data$varprop 
# PC1: 0.8527091609
# PC2: 0.1273200916
# PC3: 0.0127698676
# PC4: 0.0016227136
# PC5: 0.0006766155

# First 3 eigenfunctions
par(mfrow=c(1,3))
plot(pca_data$harmonics[1,], ylab='FPC1', ylim=c(-0.1,0.1))
plot(pca_data$harmonics[2,], ylab='FPC2', ylim=c(-0.1,0.1))
plot(pca_data$harmonics[3,], ylab='FPC3', ylim=c(-0.1,0.1))

# Interpret the PC's
par(mfrow=c(1,3))
plot(pca_data)
# 1st: Hot vs. cold areas
# 2nd: Areas where hot comes earlier vs. areas where hot comes later
# 3rd: Areas where, at the peak of the hot (and of the cold) it is hotter w.r.t. the mean
#      vs. it is colder w.r.t. the mean

# ------------------------------------------------------------------------------------
# 4. Request c)
# ------------------------------------------------------------------------------------
# Having reported a qualitative plot of the scores along the first 2 functional
# principal components, use the categorical variable zone to further enhance
# the interpretations.

data = read.table('watertemp.txt', header=T)
levels = as.factor(data$Zone)

scores_and_levels = as.data.frame(pca_data$scores[,1:2])
scores_and_levels[, 'Level'] = levels

ind_deep = which(scores_and_levels$Level == 'Deep')
ind_med = which(scores_and_levels$Level == 'Medium')
ind_surf = which(scores_and_levels$Level == 'Surface')

par(mfrow=c(1,1))
plot(scores_and_levels[,1:2])
points(scores_and_levels[ind_deep, 1:2], col='red', , pch=19)
points(scores_and_levels[ind_med, 1:2], col='blue' , pch=19)
points(scores_and_levels[ind_surf, 1:2], col='green', pch=19)
# There seem to be 3 clear clusters based on the locations.

# ------------------------------------------------------------------------------------
# 5. Request d)
# ------------------------------------------------------------------------------------
# Propose a possible dimensionality reduction for the data and discuss the results.
pca_data$varprop

# Based on the variance captured by the first PCs, we can say that we can reduce the
# dimension to 2 PCs (and it'll be enough).