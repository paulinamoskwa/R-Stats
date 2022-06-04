###########################################################################################
###                                       SPLINES                                       ###
###########################################################################################
library(splines)
load('nlr_data.rda')
data = data.frame("age"=Wage$age,"wage"=Wage$wage)    # data = [X, Y]
x    = data[,1]
y    = data[,2]

### ---------------------------------------------------------------------------------------
### General Splines
### ---------------------------------------------------------------------------------------
degree   = 15    # if = 3 -> cubic splines
model_sp = lm(y ~ bs(x, degree=degree))
x_grid   = seq(range(x)[1], range(x)[2], length.out=100)
y_pred   = predict(model_sp, list(x=x_grid), se=T)
se.bands = cbind(y_pred$fit+2*y_pred$se.fit, y_pred$fit-2*y_pred$se.fit)
plot(x, y, xlim=range(x_grid), cex=.5)
lines(x_grid, y_pred$fit, lwd=2, col='blue')
matlines(x_grid, se.bands, lwd=1, col='blue', lty=3)

#-----
# We can ad some break points (and plot them)
# br        = c(seq(min(x), mean(x), length=2), seq(mean(x)+1, max(x), 1))
# br        = c(seq(min(x), max(x), 100)
# model_cut = lm(y ~ bs(x, knots=br, degree=degree))
# ..
# knot_pred = predict(model_cut, list(x=br))
# points(br, knot_pred)
#-----

#-----
# We can either add the number of knots (where we break) or we can specify the degrees 
# of freddom knowing that: degree_of_freedom = #knots + degree. If we know where to put
# the knots then we can define the break points, otherwise we specify the number of points.
# If we want to have n knots and d degree (df = n+d) we write:
# n_knot = ..
# degree = ..
# model_cut = lm(y ~ bs(x, degree=degree, df=n_knot+degree))
# ..
# knots = attr(bs(x, degree=degree, df=n+d), 'knots')
# knots 
# knot_pred = predict(model_cut, data.frame(x=knots))
# points(br, knots_pred)
#-----

#-----
# We can specify the quantiles of the data on which we want a break instead of the knots:
# br = c(quantile(x, probs=c(0.2, 0.4, 0.6, 0.8)))
#-----

#-----
# We can specify the quantiles of the data and in addition another value
# val = ..
# br  = c(quantile(x, probs=c(0.2, 0.4, 0.6, 0.8)), val)
#-----

### ---------------------------------------------------------------------------------------
### Natural Splines
### ---------------------------------------------------------------------------------------
knots     = quantile(x, probs=c(0.1, 0.5, 0.9))
bd_knots  = quantile(x, probs=c(0.05, 0.95))
model_cut = lm(y ~ ns(x, knots=knots, Boundary.knots=bd_knots))
x_grid    = seq(range(x)[1], range(x)[2], length.out=100)
y_pred    = predict(model_cut, list(x=x_grid), se=T)
se.bands  = cbind(y_pred$fit+2*y_pred$se.fit, y_pred$fit-2*y_pred$se.fit)
plot(x, y, xlim=range(x_grid), cex=.5)
lines(x_grid, y_pred$fit, lwd=2, col='blue')
matlines(x_grid, se.bands, lwd=1, col='blue', lty=3)
knot_pred = predict(model_cut, list(x=knots))
points(knots, knot_pred, pch=19, col='red', lwd=9)
bd_knots_pred = predict(model_cut, list(x=bd_knots))
points(bd_knots, bd_knots_pred, pch=19, col='green', lwd=9)

# Note: knots and bd_knots can be values and not quantiles

### ---------------------------------------------------------------------------------------
### Smoothing Splines
### ---------------------------------------------------------------------------------------
# We can specify the degree of freedom (df) or the penalty (lambda)
fit = smooth.spline(x, y, df=61)
fit = smooth.spline(x, y, lambda=1e-10)   # higher the lambda, smoother the curve
plot(x, y, cex=.5)
lines(fit, col='blue', lwd=2)

# We can let lambda/df be choosen in an automatic way (minimizing the CV/GCV error)
fit = smooth.spline(x, y, cv=T)  # CV
fit = smooth.spline(x, y, cv=F)  # GCV
plot(x, y, cex=.5)
lines(fit, col='blue', lwd=2)

fit$lambda
fit$df

# Prediction
val = c(50.5, 60.5, 90.5)
predict(fit, val)

### ---------------------------------------------------------------------------------------
### Monotone case
### ---------------------------------------------------------------------------------------
data = subset(Orange, Orange$Tree==3)
data = data[,c(2,3)]
x    = data[,1]
y    = data[,2]

library(fda)
basis    = create.bspline.basis(range(x), breaks=x, norder=4)
cvec0    = matrix(0, basis$nbasis, 1)
Wfd0     = fd(cvec0, basis)
start_fd = fdPar(Wfd0, Lfdobj=2, lambda=1e4)
res      = smooth.monotone(x, y, start_fd)
Wfd      = res$Wfdobj
beta     = res$beta
x_grid   = seq(range(x)[1], range(x)[2], length.out=100)
y_pred   = beta[1] + beta[2]*eval.monfd(x_grid, Wfd, 0)
plot(x, y, cex=.5)
lines(x_grid, y_pred)