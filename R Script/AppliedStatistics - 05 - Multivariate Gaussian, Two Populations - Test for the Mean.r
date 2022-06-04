
# Test for Two Independent Gaussian Populations

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(car)
load("mcshapiro.test.RData")

# ------------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------------
data1 = read.table("lovingalmonds.txt", head=T)
data2 = read.table("hatingalmonds.txt", head=T)

# ------------------------------------------------------------------------------------
# 2. Assumptions and EDA 
# ------------------------------------------------------------------------------------
# Dimensions
n1 = dim(data1)[1]
n2 = dim(data2)[1]
p  = dim(data1)[2]

# Estimates
data1.mean = sapply(data1, mean)
data2.mean = sapply(data2, mean)
data1.cov  = cov(data1)
data2.cov  = cov(data2)
Sp = ((n1-1)*data1.cov + (n2-1)*data2.cov)/(n1+n2-2)
Sp.inv = solve(Sp)

# Comparison of the matrices
list(S1=data1.cov, S2=data2.cov, Spooled=Sp)

# ------------------------------------------------------------------------------------
# 3. Test for the mean of a multivariate (p>1) Gaussian (two population)
# ------------------------------------------------------------------------------------
# H0: mu1 == mu2
# H1: mu1 != mu2

# Parameters
alpha = 0.01
delta0 = c(0,0)

# Test statistic
T2 = n1*n2/(n1+n2) * (data1.mean-data2.mean-delta0) %*% Sp.inv %*%
  (data1.mean-data2.mean-delta0)

# Radius of the ellpisoid
cfr.fisher = (p*(n1+n2-2)/(n1+n2+1-p))*qf(1-alpha, p, n1+n2-1-p)

# Test
T2 < cfr.fisher

# p-value
P = 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P

# Simultaneous T2 intervals
IC.T2.p1 = c(data1.mean[1]-data2.mean[1]-sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)),
             data1.mean[1]-data2.mean[1]+sqrt(cfr.fisher*Sp[1,1]*(1/n1+1/n2)))
IC.T2.p2 = c(data1.mean[2]-data2.mean[2]-sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)),
             data1.mean[2]-data2.mean[2]+sqrt(cfr.fisher*Sp[2,2]*(1/n1+1/n2)))
IC.T2 = rbind(IC.T2.p1, IC.T2.p2)
dimnames(IC.T2)[[2]] <- c('inf','sup')
IC.T2

# Comment
# Do the interals (separately) contains the value 0?
