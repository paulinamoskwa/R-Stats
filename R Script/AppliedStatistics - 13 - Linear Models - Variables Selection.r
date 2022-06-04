
# Linear Models

# ------------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------------
library(MASS)
library(car)
library(rgl)
library(leaps)
library(ISLR)

# ------------------------------------------------------------------------------------
# Part 3 - Variables Selection
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# 1. Dataset
# ------------------------------------------------------------------------------------
data = Hitters
data = na.omit(data)
head(data)

# Re-name
n = dim(data)[1]
p = dim(data)[2]-1
v = 'X1'
if(p>1){
  for(i in 2:p){
    v = c(v, paste('X',i,sep=''))
  }
}
v = c(v, 'Y')
v[19] = 'Y'
v[20] = 'X19'
v
colnames(data) = v
head(data)
attach(data)

# ------------------------------------------------------------------------------------
# 2. Best Subset Selection (exhaustive search)
# ------------------------------------------------------------------------------------
# Best Subset Selection
regfit.full = regsubsets(Y~., data=data)
summary(regfit.full)

# Best Subset Selection: we say when we stop
regfit.full = regsubsets(Y~., data=data, nvmax=19)
summary(regfit.full)

reg.summary = summary(regfit.full)

# Which one we choose:
reg.summary$which

# R-squared
reg.summary$rsq

# R.adj^2
reg.summary$adjr2

# SSres (residual sum of squares)
reg.summary$rss

# Plots
par(mfrow=c(1,3))
plot(reg.summary$rsq, xlab="Number of Variables", ylab="R-squared", type="b")
plot(reg.summary$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="b")
plot(reg.summary$rss, xlab="Number of Variables", ylab="RSS", type="b")

# We want the model with max r.adj^2 so we extract the coefficients of that model
# Note: ind = how many coefficients has the model
ind = which.max(reg.summary$adjr2)
coef(regfit.full, ind)

# Graphical table of best results
plot(regfit.full, scale="r2", main="Exhaustive search")
plot(regfit.full, scale="adjr2", main="Exhaustive search")

# ------------------------------------------------------------------------------------
# 3. Forward and Backward Stepwise Selection
# ------------------------------------------------------------------------------------
# Forward
regfit.fwd = regsubsets(Y~.,data=data,nvmax=19,method="forward")
summary(regfit.fwd)

# Plot
par(mfrow=c(1,3))
plot(summary(regfit.fwd)$rsq, xlab="Number of Variables", ylab="R-squared", type="b")
plot(summary(regfit.fwd)$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="b")
plot(summary(regfit.fwd)$rss, xlab="Number of Variables", ylab="RSS", type="b")

par(mfrow=c(1,2))
plot(regfit.fwd,scale="r2",main="Forward Stepwise Selection")
plot(regfit.fwd,scale="adjr2",main="Forward Stepwise Selection")

# Backward
regfit.bwd = regsubsets(Y~.,data=data,nvmax=19,method="backward")
summary(regfit.bwd)

# Plot
par(mfrow=c(1,3))
plot(summary(regfit.bwd)$rsq, xlab="Number of Variables", ylab="R-squared", type="b")
plot(summary(regfit.bwd)$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="b")
plot(summary(regfit.bwd)$rss, xlab="Number of Variables", ylab="RSS", type="b")

par(mfrow=c(1,2))
plot(regfit.bwd,scale="r2",main="Backward Stepwise Selection")
plot(regfit.bwd,scale="adjr2",main="Backward Stepwise Selection")

# ------------------------------------------------------------------------------------
# 4. Comparison
# ------------------------------------------------------------------------------------
coef(regfit.full,7) # Exhaustive search
coef(regfit.fwd,7) # Forward Stepwise Selection
coef(regfit.bwd,7) # Backward Stepwise Selection

# ------------------------------------------------------------------------------------
# 5. Choosing among models using the k-fold cross-validation approach
# (exhaustive search)
# ------------------------------------------------------------------------------------
k = 10
folds = sample(1:k,nrow(data),replace=TRUE)
table(folds)

# Function that performs the prediction for regsubsets
predict.regsubsets = function(object,newdata,id){
  form  = as.formula(object$call[[2]])
  mat   = model.matrix(form,newdata)
  coefi = coef(object,id=id)
  xvars = names(coefi)
  mat[,xvars]%*%coefi
}

cv.errors = matrix(NA,k,19, dimnames=list(NULL, paste(1:19)))

for(j in 1:k){
  best.fit = regsubsets(Y~.,data=data[folds!=j,],nvmax=19)
  for(i in 1:19){
    pred = predict(best.fit,data[folds==j,],id=i)
    cv.errors[j,i] = mean( (data$Y[folds==j]-pred)^2 )
  }
}

cv.errors

root.mean.cv.errors = sqrt(apply(cv.errors,2,mean)) # average over the columns
root.mean.cv.errors

# Plot
plot(root.mean.cv.errors,type='b')
points(which.min(root.mean.cv.errors),
       root.mean.cv.errors[which.min(root.mean.cv.errors)], col='red',pch=19)
which.min(root.mean.cv.errors)

# Estimation on the full dataset
reg.best = regsubsets(Y~.,data=data, nvmax=19)
coef(reg.best,10)
