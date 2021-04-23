cureEM <- function(formula, cureformula, data, id = NULL, weights, corstr = "ind") {
  call <- match.call()
  data <- data[complete.cases(data),]
  n <- nrow(data)
  mf <- model.frame(formula, data)
  cvars <- all.vars(cureformula)
  Y <- model.extract(mf, "response")
  X <- model.matrix(attr(mf, "terms"), mf)
  Z <- model.matrix(cureformula, data)
  ## Don't need to remove intercept for ls approach
  ## X <- X[,-(colnames(X) %in% "(Intercept)")]
  ## Z <- Z[,-(colnames(Z) %in% "(Intercept)")]
  id <- eval(call["id"][[1]], data)
  if (is.null(id)) id <- 1:nrow(data)
  wgt <- eval(call["weights"][[1]], data)
  if (is.null(wgt)) wgt <- rep(1, nrow(data))
  ## initial values
  fit.a <- glm(Y[,2] ~ Z - 1, family = binomial)
  fit.b <- aftgee(Y ~ X - 1, id = id, B = 0, weights = wgt, corstr = corstr)
  emfit <- fpiter(par = c(coef(fit.a), coef(fit.b),
                          aftsurv(Y, X, coef(fit.b), wgt)),
                  X = X, Z = Z, Y = Y, id = id, weights = wgt,
                  corstr = corstr,
                  fixptfn = oneStep, 
                  control = list(maxiter = 100, tol = 1e-4))
  list(cureCoef = emfit$par[1:NCOL(Z)],
       aftCoef = emfit$par[1:NCOL(X) + NCOL(Z)],
       surv = emfit$par[-c(1:NCOL(X) + NCOL(Z))],
       iter = emfit$iter,
       convergence = emfit$convergence)
}

## EM algorithm based on the estimated regression coefficients and S(.)
oneStep <- function(para, X, Z, Y, id, weights, corstr) {
  ## extract from `para`
  a <- para[1:NCOL(Z)]
  b <- para[1:NCOL(X) + NCOL(Z)]
  sp <- para[-c(1:(NCOL(X) + NCOL(Z)))]
  pi <- 1 / (1 + c(exp(-Z %*% a)))
  w <- ifelse(Y[,2] > 0, 1, pi * sp / (1 - pi + pi * sp))
  fit.a <- glm(w ~ Z - 1, weights = weights, family = quasibinomial, control = list(maxit = 100))
  fit.b <- aftgee(Y ~ X - 1, id = id, weights = w * weights, B = 0, corstr = corstr)
  return(c(coef(fit.a), coef(fit.b), aftsurv(Y, X, coef(fit.b), w * weights)))
}

aftsurv <- function(Y, X, b, weights) {
  err <- log(Y[,1]) - X %*% b
  ## km <- survfit(Surv(err, Y[,2]) ~ 1, weights = weights)
  ## approx(x = sort(Y[,1]), y = km$surv, Y[,1])$y
  km2 <- survfit(Surv(err, Y[,2]) ~ 1, weights = weights, subset = Y[,2] > 0)
  approx(x = sort(Y[err %in% km2$time,1]), y = km2$surv, Y[,1], yleft = 1, yright = 0)$y
}