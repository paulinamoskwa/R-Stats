###########################################################################################
###                           CONFORMAL PREDICTION INTERVALS                            ###
###                             library(conformalInference)                             ###
###########################################################################################

# To install the library:
#   library(devtools)
#   install_github(repo="ryantibs/conformal", subdir="conformalInference")

library(conformalInference)
seed  = 100
B     = 1000
alpha = 0.1

# The following are full conformal. For split: conformal.pred.split()

### ---------------------------------------------------------------------------------------
### Univariate
### ---------------------------------------------------------------------------------------
load('uphillperformance.rda')
data   = df_2
y      = data$u.speed
x      = data$t.days
x_grid = seq(range(x)[1], range(x)[2], length.out =50)

### |--------|
### | Linear |
### |--------|
fit   = lm(y ~ x)
preds = predict(fit, list(x=x_grid), se=T)
plot(x, y, xlim=range(x_grid), cex=.5, pch=19)
lines(x_grid, preds$fit, lwd=2, col='blue')

# Conformal Interval evaluated in new_val
new_val    = mean(x)
lm_train   = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun
c_preds    = conformal.pred(x, y, new_val, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
cbind(c_preds$lo, c_preds$pred, c_preds$up)

# Conformal Intervals plotted
c_pr_plot = conformal.pred(x, y, x_grid, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
lines(x_grid, c_pr_plot$pred, lwd=3, col='red', lty=3)
matlines(x_grid, cbind(c_pr_plot$up, c_pr_plot$lo), lwd=3, col='red', lty=3)

### |------------|
### | Polynomial |
### |------------|
fit   = lm(y ~ poly(x, degree=2))
preds = predict(fit, list(x=x_grid), se=T)
plot(x, y, xlim=range(x_grid), cex=.5, pch=19)
lines(x_grid, preds$fit, lwd=2, col='blue')

# Conformal Interval evaluated in new_val
design_mx  = matrix(poly(x, degree=2), ncol=2)
point_0    = mean(x)
new_val    = matrix(poly(point_0, degree=2, coefs=attr(poly(x, degree=2), 'coefs')), ncol=2)
lm_train   = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun
c_preds    = conformal.pred(design_mx, y, new_val, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
cbind(c_preds$lo, c_preds$pred, c_preds$up)

# Conformal Intervals plotted
x_pred    = matrix(poly(x_grid, degree=2, coefs=attr(poly(x, degree=2), 'coefs')), ncol=2)
c_pr_plot = conformal.pred(design_mx, y, x_pred, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
lines(x_grid, c_pr_plot$pred, lwd=3, col='red', lty=3)
matlines(x_grid, cbind(c_pr_plot$up, c_pr_plot$lo), lwd=3, col='red', lty=3)

### |---------|
### | Splines |
### |---------|
library(splines)
br    = c(quantile(x, probs=c(0.25, 0.5, 0.75)))
fit   = lm(y ~ bs(x, degree=3, knots=br))
preds = predict(fit, list(x=x_grid), se=T)
plot(x, y, xlim=range(x_grid), cex=.5, pch=19)
lines(x_grid, preds$fit, lwd=2, col='blue')

# Conformal Interval evaluated in new_val
design_mx  = bs(x, degree=3, knots=br)
point_0    = mean(x)
new_val    = matrix(bs(point_0, degree=3, knots=br), nrow=1)
lm_train   = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun
c_preds    = conformal.pred(design_mx, y, new_val, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
cbind(c_preds$lo, c_preds$pred, c_preds$up)

# Conformal Intervals plotted
x_pred    = matrix(bs(x_grid, degree=3, knots=br), nrow=length(x_grid))
c_pr_plot = conformal.pred(design_mx, y, x_pred, alpha=alpha, train.fun=lm_train, predict.fun=lm_predict)
lines(x_grid, c_pr_plot$pred, lwd=3, col='red', lty=3)
matlines(x_grid, cbind(c_pr_plot$up, c_pr_plot$lo), lwd=3, col='red', lty=3)

### |-------------------|
### | Smoothing Splines |
### |-------------------|
fit = smooth.spline(x, y, cv=T)
opt = fit$df
plot(x, y, cex=.5, pch=19)
lines(fit, col='blue', lwd=2)

# Conformal Interval evaluated in new_val
new_val    = mean(x)
train_ss   = function(x,y,out=NULL){smooth.spline(x,y,df=opt)}
predict_ss = function(obj, new_x){predict(obj,new_x)$y}
c_preds    = conformal.pred(x, y, new_val, alpha=alpha, train.fun=train_ss, predict.fun=predict_ss)
cbind(c_preds$lo, c_preds$pred, c_preds$up)

# Conformal Intervals plotted
c_pr_plot = conformal.pred(x, y, x_grid, alpha=alpha, train.fun=train_ss, predict.fun=predict_ss)
lines(x_grid, c_pr_plot$pred, lwd=3, col='red', lty=3)
matlines(x_grid, cbind(c_pr_plot$up, c_pr_plot$lo), lwd=3, col='red', lty=3)

### ---------------------------------------------------------------------------------------
### Multivariate (GAM)
### ---------------------------------------------------------------------------------------
library(mgcv)
library(rgl)
load('uphillperformance.rda')
data        = df_2
y           = data$u.speed
x1          = data$t.days
x2          = data$w.ski
x1_grid     = seq(range(x1)[1], range(x1)[2], length.out =50)
x2_grid     = seq(range(x2)[1], range(x2)[2], length.out =50)
grid        = expand.grid(x1_grid, x2_grid)
names(grid) = c('x1', 'x2')

fit  = gam(y ~ s(x1, bs='cr') + s(x2, bs='cr'))
pred = predict(fit, newdata=grid)
persp3d(x1_grid, x2_grid, pred, col='yellow')
points3d(x1, x2, y, col='black', size=5)

train_gam = function(x, y, out=NULL){
  colnames(x) = c('x1', 'x2')
  train_data = data.frame(y,x)
  fit = gam(y ~ s(x1, bs='cr') + s(x2, bs='cr'), data=train_data)
}

predict_gam = function(obj, new_x){
  new_x = data.frame(new_x)
  colnames(new_x) = c('x1', 'x2')
  predict.gam(obj, new_x)
}

c_preds = conformal.pred(cbind(x1,x2), y, c(mean(x1),mean(x2)), alpha=alpha, train.fun=train_gam, predict.fun=predict_gam)
cbind(c_preds$lo, c_preds$pred, c_preds$up)