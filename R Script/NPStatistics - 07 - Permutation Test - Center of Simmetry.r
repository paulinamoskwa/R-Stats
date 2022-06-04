###########################################################################################
###                                  PERMUTATION TEST                                   ###
###                                 Center of simmetry                                  ###
###########################################################################################
seed = 100
B    = 1000

### ---------------------------------------------------------------------------------------
### One dimensional 
### ---------------------------------------------------------------------------------------
### Two-sided test
###    H0: "center" is the center of simmetry
###    H1: "center" is *not* the center of simmetry
### Test statistic: distance between the sample mean/median and the center 
data   = cars$dist
n      = length(data)
center = c(42)
boxplot(data)
plot(data, cex=0.5, pch=19)
points(center, pch=19, col='red')

# CMC to estimate the p-value 
set.seed(seed)
data_mean = mean(data)    # or: colMedianians(data)
T0        = abs(data_mean - center)
T_stat    = numeric(B)
for(k in 1:B){
  signs.perm     = rbinom(n, 1, 0.5)*2-1
  data_perm      = center + (data-center)*signs.perm
  data_perm_mean = mean(data_perm)
  T_stat[k]      = abs(data_perm_mean-center)
}

# Permutational distribution of T
hist(T_stat, xlim=range(c(T_stat,T0)), breaks=30)
abline(v=T0, col=3, lwd=2)
plot(ecdf(T_stat), xlim=range(c(T_stat,T0)))
abline(v=T0, col=3, lwd=2)

# p-value
p_value = sum(T_stat>=T0)/B
p_value

### ---------------------------------------------------------------------------------------
### Multi dimensional
### ---------------------------------------------------------------------------------------
### Two-sided test
###    H0: "center" is the center of simmetry
###    H1: "center" is *not* the center of simmetry
### Test statistic: squared distance between the sample mean/median and the center (other on notes)
data   = cars
n      = dim(data)[1]
p      = dim(data)[2]
center = c(15,1)
boxplot(data)
plot(data, cex=0.5, pch=19)
points(center[1], center[2], pch=19, col='red')

# CMC to estimate the p-value 
set.seed(seed)
data_mean = colMeans(data)    # or: colMedianians(data)
T0        = as.numeric((data_mean-center) %*% (data_mean-center))
T_stat    = numeric(B)
for(k in 1:B){
  signs.perm     = rbinom(n, 1, 0.5)*2-1
  data_perm      = center + (data-center)*matrix(signs.perm, nrow=n, ncol=p, byrow=FALSE)
  data_perm_mean = colMeans(data_perm)
  T_stat[k]      = as.numeric((data_perm_mean-center) %*% (data_perm_mean-center))
}

# Permutational distribution of T
hist(T_stat, xlim=range(c(T_stat,T0)), breaks=30)
abline(v=T0, col=3, lwd=2)
plot(ecdf(T_stat), xlim=range(c(T_stat,T0)))
abline(v=T0, col=3, lwd=2)

# p-value
p_value = sum(T_stat>=T0)/B
p_value