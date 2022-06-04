###########################################################################################
###                                   SIGNED RANK TEST                                  ###
###                              Wilcoxon Signed Rank W test                            ###
###########################################################################################
seed = 100
B    = 1000

### ---------------------------------------------------------------------------------------
### One sample
### ---------------------------------------------------------------------------------------
data    = cars$speed
q       = 0.5
c       = 10
n       = length(data)
ranks   = rank(abs(data-c))
W.plus  = sum(ranks[data>c])
W.minus = sum(ranks[data<c])   # "W.plus + W.minus" should be equal to "n*(n+1)/2
W       = W.plus - W.minus

# MC computation of the p-value
set.seed(seed)
W.sim = numeric(B)
for(k in 1:B){
  ranks.temp = sample(1:n)
  signs.temp = 2*rbinom(n, 1, q) - 1
  W.temp     = sum(signs.temp*ranks.temp)
  W.sim[k]   = W.temp
}

hist(W.sim, xlim=c(-n*(n+1)/2, n*(n+1)/2))
abline(v=W, col='red', lwd=3)
abline(v=0, lwd=3)

### Two-sided test ------------------------------------------------------------------------
###    H0: P(X>c)  = q  
###    H1: P(X>c) != q  
p_value = 2*sum(W.sim>= abs(W))/B
p_value

### One-side test -------------------------------------------------------------------------
###    H0: P(X>c) = q
###    H1: P(X>c) > q
p_value = sum(W.sim>=abs(W))/B
p_value

### One-side test -------------------------------------------------------------------------
###    H0: P(X>c) = q
###    H1: P(X>c) < q
p_value = sum(W.sim<=abs(W))/B
p_value

### ---------------------------------------------------------------------------------------
### Two samples-paired
### ---------------------------------------------------------------------------------------
data    = cars                      # = [X,Y]
diff    = data[,1] - data[,2]       # = X - Y
boxplot(diff)
q       = 0.5
n       = length(diff)
ranks   = rank(abs(diff))
W.plus  = sum(ranks[diff>0])
W.minus = sum(ranks[diff<0])   # "W.plus + W.minus" should be equal to "n*(n+1)/2
W       = W.plus - W.minus

# MC computation of the p-value
set.seed(seed)
W.sim = numeric(B)
for(k in 1:B){
  ranks.temp = sample(1:n)
  signs.temp = 2*rbinom(n, 1, q) - 1
  W.temp     = sum(signs.temp*ranks.temp)
  W.sim[k]   = W.temp
}

hist(W.sim, xlim=c(-n*(n+1)/2, n*(n+1)/2))
abline(v=W, col='red', lwd=3)
abline(v=0, lwd=3)

### Two-sided test ------------------------------------------------------------------------
###    H0: P(X>Y)  = q    ==>    H0: P(Z>0)  = q
###    H1: P(X>Y) != q           H1: P(Z>0) != q
p_value = 2*sum(W.sim>= abs(W))/B
p_value

### One-side test -------------------------------------------------------------------------
###    H0: P(X>Y) = q    ==>    H0: P(Z>0) = q
###    H1: P(X>Y) > q           H1: P(Z>0) > q
p_value = sum(W.sim>=abs(W))/B
p_value

### One-side test -------------------------------------------------------------------------
###    H0: P(X>Y) = q    ==>    H0: P(Z>0) = q    
###    H1: P(X>Y) < q           H1: P(Z>0) < q    
p_value = sum(W.sim<=abs(W))/B
p_value