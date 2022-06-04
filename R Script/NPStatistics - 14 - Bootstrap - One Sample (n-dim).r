###########################################################################################
###                                      BOOTSTRAP                                      ###
###                                 One sample (n-dim)                                  ###
###########################################################################################
seed  = 100
B     = 1000
alpha = 0.1

load('snowlevel.rda')
data = df_1[which(df_1$region=='Aosta_Valley' & df_1$town<10), c(1,2)]  # data = [feature_1, feature_2]

### |---------------------------------|
### | Bootstrap: bidimensional median |
### |---------------------------------|
library(DepthProc)
method_name = 'Tukey'

set.seed(seed)
T0     = depthMedian(data, depth_params=list(method=method_name))
T_boot = cbind(numeric(B), numeric(B))
for(k in 1:B){
  ind        = sample(1:dim(data)[1], replace=T)
  data_boot  = data[ind,]
  T_boot[k,] = depthMedian(data_boot, depth_params=list(method=method_name))
}

### |-------------------|
### | Bias and Variance |
### |-------------------|
apply(T_boot, 2, var)           # Variance
apply(T_boot, 2, mean) - T0     # Bias

### |----------------------|
### | Confidence Intervals |
### |----------------------|
# Reverse Percentile - 1
right_quantile = quantile(T_boot[,1], 1-alpha/2)
left_quantile  = quantile(T_boot[,1], alpha/2)
CI_1 = c(T0[1]-(right_quantile-T0[1]), T0[1], T0[1]-(left_quantile-T0[1]))
CI_1

# Reverse Percentile - 2
right_quantile = quantile(T_boot[,2], 1-alpha/2)
left_quantile  = quantile(T_boot[,2], alpha/2)
CI_2 = c(T0[2]-(right_quantile-T0[2]), T0[2], T0[2]-(left_quantile-T0[2]))
CI_2