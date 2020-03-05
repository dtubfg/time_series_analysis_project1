#' ---
#' title: "Exam Project 1 - 02417 Time Series Analysis"
#' author: "Christopher Scott Warhuus (s153465)"
#' date: "Due: March 7th, 2020"
#' 
#' output:
#'   pdf_document:
#'    includes:
#'      in_header: preamble-latex.tex
#' ---
#' 
## ----setup, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
library(bookdown)

#' 
#' ## Introduction
#' 
#' In this project, we attempt to fit several linear models to a time series of observations of atmospheric CO2 from Mauna Loa, Hawaii. The data is provided by the NOAA.\footnote{https://www.esrl.noaa.gov/gmd/ccgg/trends/} We include observations from January, 1970 to January, 2020 (601 observations in total).
#' 
#' ## Question 1.1: Plotting
#' 
#' A plot of the CO2 as a function of time is shown in the figure below. The training data includes all observations through 2017 and the testing data includes all observations thereafter.
#' 
## ----data, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------
df <- data.frame(read.table("/Users/Jan Warhuus/Desktop/DTU/6_semester/time_series_analysis/A1_co2.txt", header=TRUE))
N <- nrow(df)
df$new_time = c(1:N)
N_train = length(df$year[df$year <= 2017])
N_test = N - N_train
train <- df[order(df$new_time),][c(1:N_train),]
test <- df[order(df$new_time),][c((N-(N_test-1)):N),]
par(mfrow=c(1,1))
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(320,420), xlim=c(1, max(df$new_time)), lty=1)
lines(test$new_time, test$co2, col="blue")
legend(x="topleft", legend=c("Training Data", "Test Data"),
       col=c("black", "blue"), lty=1:1, cex=0.8)

#' 
#' The plot has two major components: a cyclical, or seasonal component, where CO2 levels are higest in spring/early summer and lowest in fall/early winter, and an upward trend, from about 325 ppm in 1970 to 413 ppm in 2020. The upward trend seems linear - i.e., CO2 levels are increasing at a relatively constant rate.
#' 
#' ## Question 1.2: OLS and WLS
#' 
#' We formulate a model that contains an intercept and slope - which accounts for the upward trend - as well as a cyclical component. We use a period of $p = 12$, since there are 12 observations per annual cycle (one per month). The model is:
#' 
#' $$Y_t = \alpha + \beta_t t + \beta_s \sin \left(\pi / 6\cdot t\right) + \beta_c \cos \left(\pi/6\cdot t\right) + \varepsilon_t$$
#' We estimate the parameters by solving the normal equation in theorem 3.1 (equations 3.34 and 3.35). We find the standard error of the estimated paramters using theroem 3.2, iii. The estimates can be found in the table below:
#' 
## ----least_square, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
p <- 12
x <- matrix(c(rep(1, N_train), train$new_time, sin(train$new_time * 2*pi/p), cos(train$new_time * 2*pi/p)), nrow=N_train, ncol=4)
Y <- c(train$co2)
theta_hat_ols <- solve(t(x) %*% x) %*% t(x) %*% Y

# Get variance of betas from here
# analysis <- lm(co2 ~ new_time + sin(2*pi/p*new_time) + cos(2*pi/p*new_time), data=train)
# summary(analysis)$coefficients

#Get coef info out:
errs <- (x %*% theta_hat_ols) - Y
sigma2_hat <- (t(errs) %*% errs) / (N_train - 4)
var_theta_hat_ols <- diag(sigma2_hat[1] * solve(t(x) %*% x))
tbl <- data.frame("Estimate"=theta_hat_ols, "Std deviation"=sqrt(var_theta_hat_ols))
rownames(tbl) <- c("$\\alpha$", "$\\beta_t$", "$\\beta_s$", "$\\beta_c$")
kable(tbl, caption="OLS estimates of the parameters", col.names = gsub("[.]", " ", names(tbl)))

#' 
#' Plotting the OLS model along with the data yields:
#' 
## ----OLS plot, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(320,420), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(train$new_time, x %*% theta_hat_ols, col="green")

legend(x="topleft", legend=c("Training Data", "Testing Data", "OLS Fit"),
       col=c("black", "blue", "green"), lty=1:1, cex=0.8)

#' 
#' The model fits the data relatively well, though we see slight deficiencies as the slope shifts. This particularaly pronounced since 2015, where the trend has increased slightly. We could inspect the model errors by plotting them, but we can already see that they are not independent of one another; virtually all errors in the last five years have been negative, while the errors near month 300 have been mostly positive.
#' 
#' We'll now assume the given correlation structure of the residuals to construct a non-identity $\Sigma$ matrix for the WLS model. It should be noted that, although we now consider the residuals correlated, we do not give greater weight to any particular observations in time.
#' 
#' To find $\Sigma$, we apply the relaxation algorithm (page 41). We start by finding $\rho = Cor(\varepsilon_{t_{i}},\varepsilon_{t_{j}}),$ where $|i - j| = 1$. This is equivalent to the autocorrelation of the residuals with a lag of 1, so we use R's autocorrelation function, \texttt{acf}. With this estimate of $\rho$, we construct $\Sigma$. Finally, we use $\Sigma$ to obtain the WLS estimate of the parameters (theorem 3.3, equation 3.40) and, from there, calculate the model errors and update $\rho$. We run the algorithm 5 times. The iterates of $\rho$ are given in the R output below:
#' 
## ----relaxation, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
# Take the residuals from last time
res <- errs
rho <- acf(res, lag.max=1, plot=FALSE)[[1]][2]
Sigma <- matrix(NA, ncol=N_train, nrow=N_train)
rhos <- c(NA, nrows=5)
theta_hat_wls <- theta_hat_ols

for (iter in 1:5) {
  for (i in 1:N_train) {
    for (j in 1:N_train) {
      if (i == j) { Sigma[i, j] = 1 }
      else { Sigma[i, j] <- rho[1]^abs(j - i) }
    }
  }
  print(rho)
  theta_hat_wls <- solve( t(x) %*% solve(Sigma) %*% x ) %*% t(x) %*% solve(Sigma) %*% Y
  res <- Y - (x %*% theta_hat_wls)
  rho <- acf(res, lag.max=1, plot=FALSE)[[1]][2]
}

#' 
#' The estimates of the parameters are given in the table below. The shown uncertainty measures are found using equation 3.34.
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get coef info out:
errs <- (x %*% theta_hat_wls) - Y
sigma2_hat <- (t(errs) %*% solve(Sigma) %*% errs) / (N_train - 4)
var_theta_hat_wls <- diag(sigma2_hat[1] * solve(t(x) %*% solve(Sigma) %*% x))
tbl <- data.frame("Estimate"=theta_hat_wls, "Std deviation"=sqrt(var_theta_hat_wls))
rownames(tbl) <- c("$\\alpha$", "$\\beta_t$", "$\\beta_s$", "$\\beta_c$")
kable(tbl, caption="WLS estimates of the parameters", col.names = gsub("[.]", " ", names(tbl)))

#' 
#' The WLS estimates of the parameters are notably similar to the OLS estimates. The standard deviations, however, are higher. This reflects the multiplication of $\Sigma^{-1}$ in equation 3.40 as compared to the equation in theorem 3.2, iii.
#' 
#' The similarity of the parameters is exemplified when plotting both models with the training and testing data, as seen below. We note, however, that the WLS model is very slightly steeper.
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(320,420), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(train$new_time, x %*% theta_hat_ols, col="green")
lines(train$new_time, x %*% theta_hat_wls, col="red")

legend(x="topleft", legend=c("Training Data","Testing Data", "OLS Fit", "WLS Fit"),
       col=c("black", "blue", "green", "red"), lty=1:1, cex=0.8)

#' 
#' ## Question 1.3: Local linear trend
#' 
#' We'll know fit a local linear model with a forgetting factor of 0.9. We define $L$ for this model by combining the transition matrices given for the linear trend and quadratic models on pages 53 and 54, respectively:
#' 
#' $$
#' \bm{L} =
#' \begin{bmatrix}
#' 1 & 0 & 0 & 0\\
#' 1 & 1 & 0 & 0\\
#' 0 & 0 & \cos(\pi/6) & \sin(\pi/6)\\
#' 0 & 0 & -\sin(\pi/6) & \cos(\pi/6)\\
#' \end{bmatrix}
#' $$
#' We define $\bm{f}$ as:
#' 
#' $$ 
#' \bm{f}(j) =
#' \begin{bmatrix}
#' 1\\
#' j\\
#' \sin(j\cdot\pi/6)\\
#' \cos(j\cdot\pi/6)
#' \end{bmatrix}
#' $$
#' Thus, $\bm{f}(0) = \begin{bmatrix} 1 & 0 & 0 & 1\end{bmatrix}^T$.
#' 
#' Since we skip the first ten observations as transient, we start by finding $\bm{F}_{10}$, $\bm{h}_{10}$ and $\bm{\hat{\theta}}_{10}$ using equation 3.99. Then, for each step from $t=11$ onward, we find $\hat{Y}_N$ and $\bm{F}_{N+1}$ and $\bm{h}_{N+1}$ using theroem 3.13, equations 3.103 and 3.104, respectively. The one-step prediction errors as well as the estimates of sigma for each observation are found in the plots below. We note that the sigma estimates were found using the global estimator in the course slides.  
#' 
## ----loclinearmodel, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------
# Define the period
period <- 12

# Define L matrix
L <- matrix(c(1,1,0,0,
              0,1,0,0,
              0,0,cos(2*pi/period),-sin(2*pi/period),
              0,0,sin(2*pi/period),cos(2*pi/period)),
              nrow=4, ncol=4)

# Define function f
f <- function(j) {
  return(c(1,
           j,
           sin(2*pi/period*j),
           cos(2*pi/period*j)))
}

# Insert knowns (p is numer of params)
Y <- c(df$co2)
time <- c(1:N)
lambda <- 0.9
first <- 10
alpha <- 0.05
p <- 4

# Effective number of observations is 10 with lambda = 0.9
T_ <- 10

# Start with 0 intial matrices
# FN <- matrix(rep(0, length(f(0))^2), nrow=length(f(0)), ncol=length(f(0)))
# hN <- rep(0, length(f(0)))

# Find Sigma
SigmaInv_generator <- function(lamb, start, n, time) {
  powers <- sort(time[start:n], decreasing = TRUE) - 1
  SigInv <- diag(x = lambda^powers, nrow=n, ncol=n)
  return(SigInv)
}

# Find x
x_generator <- function(funct, n, params) {
  x <- matrix(c(NA, n*params), nrow=n, ncol=params)
  for (i in 1:n) {
    x[i,] <- t(funct(-n + i))
  }
  return(x)
}

# Set empty vectors to save data
pred <- rep(0, N)
err <- rep(0, N)
sigma2.global <- rep(0, N)
sigma2.local <- rep(0, N)
sigma2.conf <- rep(0, N)
local_means <- rep(0, N)
local_means.conf <- rep(0, N)

# Find FN and hN for first 10 elements (which we do not find theta_hat or predictions for)
x <- x_generator(funct=f, n=first, params=p)
SigmaInverse <- SigmaInv_generator(lamb=lambda, start=1, n=first, time=time)
FN <- t(x) %*% SigmaInverse %*% x
hN <- t(x) %*% SigmaInverse %*% Y[1:first]
theta_hat_loclin <- solve(FN) %*% hN

# Keep FNs and hNs
FNs <- list(0, 0, 0, 0, 0, 0, 0, 0, 0, FN)
hNs <- list(0, 0, 0, 0, 0, 0, 0, 0, 0, hN)

# And now start performing updates
for (tN in first:(N_train-1)) {
  
  # Now make prediction and find the error on the prediction
  pred[tN + 1] <- t(f(1)) %*% theta_hat_loclin
  err[tN + 1] <- (pred[tN + 1] - Y[tN + 1])
  
  # Global error
  sig2 <- 0
  for (j in (first + 1):(tN + 1)) {
    first_term <- 1/((tN + 1) - first)
    numerator <- err[j]^2
    denominator <- 1 + (t(f(1)) %*% solve(FNs[[j-1]]) %*% f(1))
    sig2 <- sig2 + first_term * (numerator/denominator) #<- sigma estimator from slides
  } 
  sigma2.global[tN + 1] <- sig2
  
  # Confidence interval
  if (tN > first) {
    var <- sigma2.global[tN] * (1 + t(f(1)) %*% solve(FNs[[tN]]) %*% f(1))
    sigma2.conf[tN + 1] <- qt(alpha/2, T_-p) * sqrt(var) #<- equation 3.92
  }
  
  
  FNs[[tN + 1]] <- FNs[[tN]] + lambda^tN * f(-tN) %*% t(f(-tN))
  hNs[[tN + 1]] <- lambda * solve(L) %*% hNs[[tN]] + f(0) * Y[tN + 1]
  theta_hat_loclin <- solve(FNs[[tN + 1]]) %*% hNs[[tN + 1]]
  local_means[tN + 1] <- theta_hat_loclin[[1]]
  
  x <- x_generator(funct=f, n=(tN+1), params=p)
  SigmaInverse <- SigmaInv_generator(lamb=lambda, start=1, n=(tN+1), time=time)
  
  var_theta0 <- sigma2.global[tN + 1] * solve(t(x) %*% SigmaInverse %*% x)
  local_means.conf[tN + 1] <- qt(1-alpha/2, T_-p) * sqrt(var_theta0[1,1])
}

# So that was the training data... Now we want to stop updating...
l <- 1
s2 <- sigma2.global[N_train]
for (tN in N_train:(N-1)) {
  
  # Now make prediction and find the error on the prediction
  pred[tN + 1] <- t(f(l)) %*% theta_hat_loclin
  err[tN + 1] <- (pred[tN + 1] - Y[tN + 1])
  
  # Confidence interval
  var <- s2 * (1 + t(f(l)) %*% solve(FNs[[N_train]]) %*% f(l))
  sigma2.conf[tN + 1] <- qt(alpha/2, T_-p) * sqrt(var)
  l = l+1
}

par(mfrow=c(2,1), cex=0.8, mai=c(0.6,0.8,0.4,0.2))
plot(time[(first+1):N_train], sqrt(sigma2.global[(first+1):N_train]), xlab="t", ylab="sigma")
plot(time[(first+1):N_train], err[(first+1):N_train], type="l", ylab="error", xlab="t")


#' 
#' It seems the estimate for sigma converges to approximately 0.6, then follows an slight upward trend. The errors seem relatively random, though there is a bit of periodicity that our model does not pick up on.
#' 
#' The observations and local linear model estimates are plotted below. The confidence intervals were found using equation 3.92 and the estimates of sigma plotted above.
#' 
## ----loclinplot1, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(min(df$co2), max(pred+sigma2.conf)), xlim=c(1, N), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train], col="orange")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train] + sigma2.conf[(first+1):N_train], col="orange", lty="dashed")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train] - sigma2.conf[(first+1):N_train], col="orange", lty="dashed")
test_start <- N - (N_test - 1)
lines(test$new_time, pred[test_start:N], col="cyan")
lines(test$new_time, pred[test_start:N] + sigma2.conf[test_start:N], col="cyan", lty="dashed")
lines(test$new_time, pred[test_start:N] - sigma2.conf[test_start:N], col="cyan", lty="dashed")

legend(x="topleft", legend=c("Training data", "Testing data", "1-step predictions on training data", "95% PI on predictions on training data", "1-step predictions on testing data", "95% PI on predictions on testing data"), col=c("black", "blue", "orange", "orange", "cyan", "cyan"), lty=c("solid", "solid", "solid", "dashed", "solid", "dashed"), cex=0.8)

#' 
#' The model seems to work quite well. As expected, the model follows the observed data quite well with $\lambda = 0.9$ and with a lag of one time step. The prediction intervals also seem correct - the predicted data is generally between the two dotted boundaries, though it oversteps occanssionally. 
#' 
## ----loclinplot2, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
N_since2010 <- length(df$year[df$year >= 2010])
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(382, 417), xlim=c(N-N_since2010+1, N), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(train$new_time[10:N_train], pred[10:N_train], col="orange")
lines(train$new_time[10:N_train], pred[10:N_train] + sigma2.conf[10:N_train], col="orange", lty="dashed")
lines(train$new_time[10:N_train], pred[10:N_train] - sigma2.conf[10:N_train], col="orange", lty="dashed")
test_start <- N - (N_test - 1)
lines(test$new_time, pred[test_start:N], col="cyan")
lines(test$new_time, pred[test_start:N] + sigma2.conf[test_start:N], col="cyan", lty="dashed")
lines(test$new_time, pred[test_start:N] - sigma2.conf[test_start:N], col="cyan", lty="dashed")

legend(x="topleft", legend=c("Training data", "Testing data", "1-step predictions on training data", "95% PI on predictions on training data", "1-step predictions on testing data", "95% PI on predictions on testing data"), col=c("black", "blue", "orange", "orange", "cyan", "cyan"), lty=c("solid", "solid", "solid", "dashed", "solid", "dashed"), cex=0.8)

#' 
#' Well now use the training data to predict 1, 2, 6, 12 and 24 time steps ahead. To do so, we use $\bm{\hat{\theta}}_{576}$ in equation 3.101. The prediction intervals are found using equation 3.92.
## ----try, echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------

lsteps <- c(1, 2, 6, 12, 24)
steps <- length(lsteps)
new_preds <- rep(0, steps)
new_errs <- rep(0, steps)
new_confs <- rep(0, steps)
vars <- rep(0, steps)

i <- 1
for (l in lsteps) {
  new_preds[i] <- t(f(l)) %*% theta_hat_loclin #<- equation 3.101
  new_errs[i] <- (new_preds[i] - Y[N_train + l])
  
  vars[i] <- sigma2.global[N_train] * (1 + t(f(l)) %*% solve(FNs[[N_train]]) %*% f(l)) #<- equation 3.102
  new_confs[i] <- qt(alpha/2, T_-p) * sqrt(vars[i]) #<- equation 3.92
  i = i+1
}

tbl <- data.frame("Months ahead"=lsteps, "Prediction"=new_preds, "2.5 pct. PI"=new_preds-new_confs, "97.5 pct. PI"=new_preds+new_confs, "Est. sigma"=sqrt(vars), "Actual value"=Y[c(N_train + lsteps)])
kable(tbl, caption="Predictions with the local linear model for 1, 2, 6, 12 and 24 months ahead.", col.names = gsub("[.]", " ", names(tbl)))


#' 
#' The predictions are not entirely far off from the test observations, though we do see that the confidence in the prediction decreases the further we predict. This is expected, since we variance of the error increases for higher prediction steps (see equation 3.92).
#' 
#' Finally, we plot the data and the estimated mean for each time step (see the plot below). The confidence interval is quite small - i.e., the estimates of the intersection ($\alpha$) are likely quite accurate.
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot estimated mean moving mean model with confidence intervals.
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(min(df$co2), max(415)), xlim=c(1, N), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(df$new_time[(first + 1):N_train], local_means[(first + 1):N_train], type="l", col="blue")
lines(df$new_time[(first + 1):N_train], local_means[(first + 1):N_train] + local_means.conf[(first + 1):N_train], lty="dashed", col="brown")
lines(df$new_time[(first + 1):N_train], local_means[(first + 1):N_train] - local_means.conf[(first + 1):N_train], lty="dashed", col="brown")
legend(x="topleft", legend=c("Training data", "Testing data", "Local mean", "95% CI on local mean"), col=c("black", "blue", "brown", "brown"), lty=c("solid", "solid", "solid", "dashed"), cex=0.8)

#' 
#' ## Question 1.4: Optimal $\lambda$
#' 
#' To find the optimal forgetting factor, we plot the squared-errors for the local linear model as a function of lambda values:
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
Y <- c(df$co2)
lambda_iters <- seq(from = 0.8, to = .99, by = 0.01)
sum_pred_errs <- rep(0, length(lambda_iters))
count <- 1
for (lamb in lambda_iters) {
  lambda <- lamb
  skip <- 100
  alpha <- 0.05
  T_ <- 10
  
  FN <- matrix(rep(0, length(f(0))^2), nrow=length(f(0)), ncol=length(f(0)))
  hN <- rep(0, length(f(0)))
  
  for (j in 0:(skip - 1)) {
    FN = FN + (lambda^j * f(-j) %*% t(f(-j)))
    hN = hN + (lambda^j * f(-j) * Y[(skip-j)])
  }
  
  pred <- rep(0, N)
  err <- rep(0, N)
  sigma2 <- rep(0, N)
  sigma2.conf <- rep(0, N)
  for (tN in skip:(N - 1)) {
    
    theta_hat_loclin <- solve(FN) %*% hN
    pred[tN + 1] <- t(f(1)) %*% theta_hat_loclin
    err[tN + 1] <- pred[tN + 1] - Y[tN + 1]
    
    sig2 <- 0
    for (j in (skip + 1):(tN + 1)) {
      sig2 <- sig2 + 1/((tN + 1) - skip) * err[j]^2/(1 + (t(f(1)) %*% solve(FN) %*% f(1)))
    }
    sigma2[tN + 1] <- sig2
    sigma2.conf[tN + 1] <- qt(1-alpha/2, T_-length(f(0))) * sqrt(sig2)
    
    FN <- FN + lambda^tN * f(-tN) %*% t(f(-tN))
    hN <- lambda * solve(L) %*% hN + f(0) * Y[tN + 1]
  }
  sum_pred_errs[count] <- sum(t(err[(skip+1):N]) %*% err[(skip+1):N])
  count = count + 1
}
plot(lambda_iters, sum_pred_errs, ylab="SSE", xlab="lambda")

#' 
#' We see that the sum of squared-errors is lowest somewhere around 0.92 to 0.95. We repeat the above with a smaller interval for the forgetting factor:
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
Y <- c(df$co2)
lambda_iters <- seq(from = 0.92, to = .95, by = 0.002)
sum_pred_errs <- rep(0, length(lambda_iters))
count <- 1
for (lamb in lambda_iters) {
  lambda <- lamb
  skip <- 100
  alpha <- 0.05
  T_ <- 10
  
  FN <- matrix(rep(0, length(f(0))^2), nrow=length(f(0)), ncol=length(f(0)))
  hN <- rep(0, length(f(0)))
  
  for (j in 0:(skip - 1)) {
    FN = FN + (lambda^j * f(-j) %*% t(f(-j)))
    hN = hN + (lambda^j * f(-j) * Y[(skip-j)])
  }
  
  pred <- rep(0, N)
  err <- rep(0, N)
  sigma2 <- rep(0, N)
  sigma2.conf <- rep(0, N)
  for (tN in skip:(N - 1)) {
    
    theta_hat_loclin <- solve(FN) %*% hN
    pred[tN + 1] <- t(f(1)) %*% theta_hat_loclin
    err[tN + 1] <- pred[tN + 1] - Y[tN + 1]
    
    sig2 <- 0
    for (j in (skip + 1):(tN + 1)) {
      sig2 <- sig2 + 1/((tN + 1) - skip) * err[j]^2/(1 + (t(f(1)) %*% solve(FN) %*% f(1)))
    }
    sigma2[tN + 1] <- sig2
    sigma2.conf[tN + 1] <- qt(1-alpha/2, T_-length(f(0))) * sqrt(sig2)
    
    FN <- FN + lambda^tN * f(-tN) %*% t(f(-tN))
    hN <- lambda * solve(L) %*% hN + f(0) * Y[tN + 1]
  }
  sum_pred_errs[count] <- sum(t(err[(skip+1):N]) %*% err[(skip+1):N])
  count = count + 1
}
plot(lambda_iters, sum_pred_errs, xlim=c(0.92, 0.95), ylim=c(300,310), ylab="SSE", xlab="lambda")

#' 
#' Zooming in on the above plot yields:
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(lambda_iters, sum_pred_errs, xlim=c(0.93, 0.94), ylim=c(303,304), ylab="SSE", xlab="lambda")

#' 
#' It would seem that the optimal forgetting factor is around 0.936. We'll assume this is close enough and move on to find the one-step predictions with this forgetting factor. The 1-step predictions are shown in the plot below:
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
lambda <- 0.936
T_ <- 16

# Set empty vectors to save data
pred <- rep(0, N)
err <- rep(0, N)
sigma2.global <- rep(0, N)
sigma2.local <- rep(0, N)
sigma2.conf <- rep(0, N)
local_means <- rep(0, N)
local_means.conf <- rep(0, N)

# Find FN and hN for first 10 elements (which we do not find theta_hat or predictions for)
x <- x_generator(funct=f, n=first, params=p)
SigmaInverse <- SigmaInv_generator(lamb=lambda, start=1, n=first, time=time)
FN <- t(x) %*% SigmaInverse %*% x
hN <- t(x) %*% SigmaInverse %*% Y[1:first]
theta_hat_loclin <- solve(FN) %*% hN

# Keep FNs and hNs
FNs <- list(0, 0, 0, 0, 0, 0, 0, 0, 0, FN)
hNs <- list(0, 0, 0, 0, 0, 0, 0, 0, 0, hN)

# And now start performing updates
for (tN in first:(N_train-1)) {
  
  # Now make prediction and find the error on the prediction
  pred[tN + 1] <- t(f(1)) %*% theta_hat_loclin
  err[tN + 1] <- (pred[tN + 1] - Y[tN + 1])
  
  # Global error
  sig2 <- 0
  for (j in (first + 1):(tN + 1)) {
    first_term <- 1/((tN + 1) - first)
    numerator <- err[j]^2
    denominator <- 1 + (t(f(1)) %*% solve(FNs[[j-1]]) %*% f(1))
    sig2 <- sig2 + first_term * (numerator/denominator) #<- sigma estimator from slides
  } 
  sigma2.global[tN + 1] <- sig2
  
  # Confidence interval
  if (tN > first) {
    var <- sigma2.global[tN] * (1 + t(f(1)) %*% solve(FNs[[tN]]) %*% f(1))
    sigma2.conf[tN + 1] <- qt(alpha/2, T_-p) * sqrt(var) #<- equation 3.92
  }
  
  
  FNs[[tN + 1]] <- FNs[[tN]] + lambda^tN * f(-tN) %*% t(f(-tN))
  hNs[[tN + 1]] <- lambda * solve(L) %*% hNs[[tN]] + f(0) * Y[tN + 1]
  theta_hat_loclin <- solve(FNs[[tN + 1]]) %*% hNs[[tN + 1]]
  local_means[tN + 1] <- theta_hat_loclin[[1]]
  
  x <- x_generator(funct=f, n=(tN+1), params=p)
  SigmaInverse <- SigmaInv_generator(lamb=lambda, start=1, n=(tN+1), time=time)
  
  var_theta0 <- sigma2.global[[1]] * solve(t(x) %*% SigmaInverse %*% x)
  local_means.conf[tN + 1] <- qt(1-alpha/2, T_-p) * sqrt(var_theta0[1,1])
}

# So that was the training data... Now we want to stop updating...
l <- 1
s2 <- sigma2.global[N_train]
for (tN in N_train:(N-1)) {
  
  # Now make prediction and find the error on the prediction
  pred[tN + 1] <- t(f(l)) %*% theta_hat_loclin
  err[tN + 1] <- (pred[tN + 1] - Y[tN + 1])
  
  # Confidence interval
  var <- s2 * (1 + t(f(l)) %*% solve(FNs[[N_train]]) %*% f(l))
  sigma2.conf[tN + 1] <- qt(alpha/2, T_-p) * sqrt(var)
  l = l+1
}

#' 
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(train$new_time, train$co2, type="l", xlab='Time [months]', ylab='CO2 [ppm]', ylim=c(382, 417), xlim=c(N-N_since2010+1, N), lty=1)
lines(test$new_time, test$co2, col="blue")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train], col="orange")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train] + sigma2.conf[(first+1):N_train], col="orange", lty="dashed")
lines(train$new_time[(first+1):N_train], pred[(first+1):N_train] - sigma2.conf[(first+1):N_train], col="orange", lty="dashed")
test_start <- N - N_test + 1
lines(test$new_time, pred[test_start:N], col="cyan")
lines(test$new_time, pred[test_start:N] + sigma2.conf[test_start:N], col="cyan", lty="dashed")
lines(test$new_time, pred[test_start:N] - sigma2.conf[test_start:N], col="cyan", lty="dashed")

legend(x="topleft", legend=c("Training data", "Testing data", "1-step predictions on training data", "95% PI on predictions on training data", "1-step predictions on testing data", "95% PI on predictions on testing data"), col=c("black", "blue", "orange", "orange", "cyan", "cyan"), lty=c("solid", "solid", "solid", "dashed", "solid", "dashed"), cex=0.8)

#' 
#' A table of predictions for 1, 2, 6, 12 and 24 months ahead of the training data is shown below. The main difference is that the estimates for sigma are quite a bit lower than for the previous result. It seems the model performs quite decently - the actual observation falls within the prediction interval for every prediction.
#' 
## ---- echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------

lsteps <- c(1, 2, 6, 12, 24)
steps <- length(lsteps)
new_preds <- rep(0, steps)
new_errs <- rep(0, steps)
new_confs <- rep(0, steps)
vars <- rep(0, steps)

i <- 1
for (l in lsteps) {
  new_preds[i] <- t(f(l)) %*% theta_hat_loclin #<- equation 3.101
  new_errs[i] <- (new_preds[i] - Y[N_train + l])
  
  vars[i] <- sigma2.global[N_train] * (1 + t(f(l)) %*% solve(FNs[[N_train]]) %*% f(l)) #<- equation 3.102
  new_confs[i] <- qt(alpha/2, T_-p) * sqrt(vars[i]) #<- equation 3.92
  i = i+1
}

tbl <- data.frame("Months ahead"=lsteps, "Prediction"=new_preds, "2.5 pct. PI"=new_preds-new_confs, "97.5 pct. PI"=new_preds+new_confs, "Est. sigma"=sqrt(vars), "Actual value"=Y[c(N_train + lsteps)])
kable(tbl, caption="Predictions with the local linear model for 1, 2, 6, 12 and 24 months ahead.", col.names = gsub("[.]", " ", names(tbl)))


#' 
#' ## Question 1.5: Overall
#' 
#' We've now trained and tested several models on the data. Though the trend in the data at first glance seemed relatively linear, we found that the OLS and WLS models were not able to pick up on important nuances, particularly the steepening trend over the past 50-60 months. Particularly for the OLS model, the assumption of independent residuals did not seem to hold.
#' 
#' Next, we tried a local linear model with a forgetting factor of 0.9, giving the model an effective 10 observations. This model seemed to perform well for short-term predictions of up to 24 months. It was found, however, that the optimal forgetting factor was slightly higher, at about 0.936. This should be the preferred model, and we see that observations lie within the prediction interval up to 24 months ahead.
#' 
#' To improve the OLS and WLS models, it may make sense to add a quadratic component to pick up on the increase in slope for the last several observations, particularly since there's no reason to believe that the slope will decrease anytime soon. Additionally, since the errors for the local model model exhbibit uncaptured oscillations, it would make sense to take a closer look at the oscillation period and include this in the model. This counts for both the local linear and global trend models.
#' 
