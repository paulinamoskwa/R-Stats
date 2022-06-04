###########################################################################################
###                                  PERMUTATION TEST                                   ###
###                             Two independent samples (n-dim)                         ###
###########################################################################################
seed = 100
B    = 1000

### Two-sided test 
###    H0: X1  =^d X2
###    H1: X1 !=^d X2
### Test statistic: squared distance between the two sample mean/median vectors
data_1    = cars
data_2    = cars+1
n1        = dim(data_1)[1]
n2        = dim(data_2)[1]
n         = n1+n2
data_pool = rbind(data_1, data_2)

# CMC to estimate the p-value 
set.seed(seed)
data_1_mean = colMeans(data_1)    # or: colMedianians(data_1)
data_2_mean = colMeans(data_2) 
T0          = as.numeric((data_1_mean-data_2_mean) %*% (data_1_mean-data_2_mean))
T_stat      = numeric(B)
for(k in 1:B){
  perm             = sample(1:n)
  data_pool_perm   = data_pool[perm,]
  data_1_perm      = data_pool_perm[1:n1,]
  data_2_perm      = data_pool_perm[(n1+1):n,]
  data_1_perm_mean = colMeans(data_1_perm)
  data_2_perm_mean = colMeans(data_2_perm)
  T_stat[k]        = as.numeric((data_1_perm_mean-data_2_perm_mean) %*% 
                                  (data_1_perm_mean-data_2_perm_mean))
}

# Permutational distribution of T
hist(T_stat, xlim=range(c(T_stat,T0)), breaks=30)
abline(v=T0, col=3, lwd=2)
plot(ecdf(T_stat), xlim=range(c(T_stat,T0)))
abline(v=T0, col=3, lwd=2)

# p-value
p_value = sum(T_stat>=T0)/B
p_value