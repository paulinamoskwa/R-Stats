###########################################################################################
###                                      BOOTSTRAP                                      ###
###                                  Test and p-values                                  ###
###########################################################################################
seed  = 100
B     = 1000
alpha = 0.1

### One dimensional
###    H0: median(data)=0 vs. H1: median(data)!=0 
data = stabledist::rstable(1000,1.8,0)
hist(data)
boxplot(data)

library(pbapply)
library(parallel)
set.seed(seed)
T0      = median(data)
T_boot  = numeric(B)
cl      = makeCluster(2)
wrapper = function(dummy){T_boot = median(sample(data, replace=T))}
clusterExport(cl=cl, list('data'))
T_boot  = pbsapply(T_boot, wrapper, cl=cl)

# Confidence Interval
right_quantile = quantile(T_boot, 1-alpha/2)
left_quantile  = quantile(T_boot, alpha/2)
CI = c(T0-(right_quantile-T0), T0, T0-(left_quantile-T0))
CI
plot(ecdf(T_boot))
abline(v=CI, col='red')

# P-value
alpha_grid = seq(0.001, 0.5, length.out=100)
CI_calc = function(alpha){
  right_quantile = quantile(T_boot, 1-alpha/2)
  left_quantile  = quantile(T_boot, alpha/2)
  CI = c(T0-(right_quantile-T0), T0-(left_quantile-T0))
  names(CI) = c('lwr','upr')
  return(CI)
}

CI_list = pblapply(alpha_grid, CI_calc)
CI_mat  = dplyr::bind_rows(CI_list)
check   = CI_mat[,1]>0 | CI_mat[,2]<0
(alpha_grid[check])[1]   # <- p-value 