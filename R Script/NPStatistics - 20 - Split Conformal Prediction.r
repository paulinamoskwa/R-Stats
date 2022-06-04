###########################################################################################
###                              SPLIT CONFORMAL PREDICTION                             ###
###########################################################################################
seed  = 100
B     = 1000
alpha = 0.1

### ---------------------------------------------------------------------------------------
### Univariate
### ---------------------------------------------------------------------------------------
parz       = read.table('parziali.txt')
data       = parz$PI
n          = length(data)
train_prop = 0.5

### |-------------------|
### | Discrepancy-based |
### |-------------------|
set.seed(seed)
train_id   = sample(1:n, ceiling(n*train_prop), replace=F)
train_set  = data[train_id]
x          = data[-train_id]
x_grid     = seq(min(x)-0.5*diff(range(x)), max(x)+0.5*diff(range(x)), length.out=100)
p_value    = numeric(length(x_grid))
train_mean = mean(train_set)   # or median

NC = function(z_aug, i){
  abs(z_aug[i] - train_mean)
}

for(k in 1:length(x_grid)){
  x_aug  = c(x, x_grid[k])
  scores = numeric(length(x_aug))
  for(i in 1:length(scores)){
    scores[i] = NC(x_aug, i)
  }
  p_value[k] = sum(scores>=scores[length(x_aug)])/(length(x_aug))
}

# Plot of the p-values
plot(x_grid, p_value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red')
points(x, numeric(length(x)), pch=3)

# Prediction Interval
PI_grid = x_grid[which(p_value>=alpha)]
PI      = c(min(PI_grid), max(PI_grid))
PI
abline(v=PI, col='red')
points(PI_grid, numeric(length(PI_grid)), pch=16, col='red')

### ---------------------------------------------------------------------------------------
### Multivariate
### ---------------------------------------------------------------------------------------
data       = read.table('parziali.txt')
n          = dim(data)[1]
train_prop = 0.5

set.seed(seed)
train_id   = sample(1:n, ceiling(n*train_prop), replace=F)
train_set  = data[train_id,]
calib_set  = data[-train_id,]
x          = calib_set[,1]
y          = calib_set[,2]
x_grid     = seq(min(x)-0.5*diff(range(x)), max(x)+0.5*diff(range(x)), length.out=20)
y_grid     = seq(min(y)-0.5*diff(range(y)), max(y)+0.5*diff(range(y)), length.out=20)
p_value    = matrix(nrow=length(x_grid), ncol=length(y_grid))

### |---------------------------|
### | Predict a new (x,y) point |
### |---------------------------|
train_mean = colMeans(train_set)

NC = function(z_aug, i){
  sum((z_aug[i,]-train_mean)^2) # Euclidean distance
}

for(k in 1:length(x_grid)){
  for(h in 1:length(y_grid)){
    data_aug = rbind(calib_set, c(x_grid[k], y_grid[h]))
    scores   = numeric(dim(data_aug)[1])
    for(i in 1:length(scores)){
      scores[i] = NC(data_aug, i)
    }
    p_value[k,h] = sum(scores>=scores[dim(data_aug)[1]])/(dim(data_aug)[1])
  }
}

# Plot of the p-value
image(x_grid, y_grid, p_value, zlim=c(0,1), asp=1)
points(calib_set, pch=16)
contour(x_grid, y_grid, p_value, levels=alpha, add=T)
points(train_mean[1], train_mean[2], pch=4, cex=4)

### |------------|
### | Regression |
### |------------|
# Model on training (can be any kind of model, check how to extract the prediction!)
fm = smooth.spline(train_set[,1], train_set[,2], lambda=0.05)

NC = function(z_aug, i){
  abs(z_aug[i,2] - predict(fm, x = z_aug[i,1])$y)
}

for(k in 1:length(x_grid)){
  for(h in 1:length(y_grid)){
    data_aug = rbind(calib_set, c(x_grid[k], y_grid[h]))
    scores   = numeric(dim(data_aug)[1])
    for(i in 1:length(scores)){
      scores[i] = NC(data_aug, i)
    }
    p_value[k,h] = sum(scores>=scores[dim(data_aug)[1]])/(dim(data_aug)[1])
  }
}

# Plot of the p-value
image(x_grid, y_grid, p_value, zlim=c(0,1))
points(data, pch=16)
contour(x_grid, y_grid, p_value, levels=alpha, add=T)