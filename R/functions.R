#
# This is the implementation in the original glmperm package (GPL-3)
# https://cran.r-project.org/web/packages/glmperm/index.html
#

summary.prr.test <-
  function(object, digits = max(3, getOption("digits") - 3),...)
  {
    cat("\n\n    Permutation of Regressor Residuals Test:\n\n\n")
    cat("Call: \n", deparse(object$call), "\n\n")
    cat("number of observations used: ",  format(round(object$nobs, digits)),"\n\n")
    cat("null hypothesis: regression coefficient of covariate", object$var,"= 0")
    cat("\nobserved Likelihood Ratio Test Statistics: ", format(round(object$LRstat, digits)),"\n\n")
    cat("    -----------------------------------------\n")
    cat("    Results based on chi-squared distribution \n")
    cat("    -----------------------------------------\n")
    cat("\nobserved p-value:", format(round(object$p.value.obs, digits)),"\n\n")
    cat("    ---------------------------------------------------\n")
    cat("    Results based on permutation of regressor residuals \n")
    cat("    ---------------------------------------------------\n")
    cat("\npermutation p-value for simulated p-values <= observed p-value:",  format(round(object$p.value.perm$p0, digits)),"(Std.err:", paste(format(round(object$p.value.perm.se$se.p0,digits)),")",sep=""),"\n")
    cat("\npermutation p-value for simulated p-values <= 1.005 observed p-value:",  format(round(object$p.value.perm$p005, digits)),"(Std.err:", paste(format(round(object$p.value.perm.se$se.p005,digits)),")",sep=""),"\n")
    cat("\npermutation p-value for simulated p-values <= 1.01 observed p-value:",  format(round(object$p.value.perm$p01, digits)),"(Std.err:", paste(format(round(object$p.value.perm.se$se.p01,digits)),")",sep=""),"\n")
    cat("\npermutation p-value for simulated p-values <= 1.02 observed p-value:",  format(round(object$p.value.perm$p02, digits)),"(Std.err:", paste(format(round(object$p.value.perm.se$se.p02,digits)),")",sep=""),"\n")
    cat("\npermutation p-value for simulated p-values <= 1.04 observed p-value:",  format(round(object$p.value.perm$p04, digits)),"(Std.err:", paste(format(round(object$p.value.perm.se$se.p04,digits)),")",sep=""),"\n\n")
    ## Warnings for Over and Underdispersion
    if(object$fit1$family$family=="binomial" & object$estimated.Dispersion >1.5){
      cat("*************************************************************************\n")
      cat("WARNING: estimated dispersion is > 1.5, rather use family = quasibinomial\n")
      cat("*************************************************************************\n")}
    if(object$fit1$family$family=="binomial" & object$estimated.Dispersion<0.5){
      cat("*************************************************************************\n")
      cat("WARNING: estimated dispersion is < 0.5, rather use family = quasibinomial\n")
      cat("*************************************************************************\n")}
    if(object$fit1$family$family=="poisson" & object$estimated.Dispersion>1.5){
      cat("************************************************************************\n")
      cat("WARNING: estimated dispersion is > 1.5, rather use family = quasipoisson\n")
      cat("************************************************************************\n")}
    if(object$fit1$family$family=="poisson" & object$estimated.Dispersion<0.5){
      cat("************************************************************************\n")
      cat("WARNING: estimated dispersion is < 0.5, rather use family = quasipoisson\n")
      cat("************************************************************************\n")}
    invisible(object)
  }


Cox <- function() {
  structure(list(family = "Cox"), class = "family")
}

glm.perm <-
  function(y, x, Family)
  {
    ## resample the data
    sel <- sample(1:nrow(x), replace = FALSE)
    x.resample <- cbind(x[, -which(colnames(x)=="resid"), drop = FALSE], x[, "resid", drop = FALSE][sel])
    ## calculate deviances and dispersion factors for resampled data
    if(Family$family == "Cox"){
      f1 <- try(coxph.fit(x=x.resample, y=y, strata=NULL, control=coxph.control(), method="efron", rownames=NULL))
      devi <- -2*f1$loglik[2]
      return(devi)
    } else {
      f1 <- glm.fit(x.resample, y, family=Family)
      df.r <- f1$df.residual
      if(f1$family$family %in% c("poisson", "binomial")){disp <- 1}
      if(all(f1$family$family != c("poisson", "binomial"))){
        if(df.r > 0){
          if(any(f1$weights == 0)){warning("observations with zero weight not used for calculating dispersion")}
          disp <- sum((f1$weights * f1$residuals^2)[f1$weights > 0])/df.r
        }
        if(df.r==0){disp <- NaN}
      }
      devi <- f1$deviance
      return(c(devi, disp))
    }
  }


prr.test <-
  function(formula, var, family=gaussian, data, nrep = 1000, seed=12345, Silent=TRUE, weights,  subset, na.action,
           start = NULL, etastart, mustart, offset, control = glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)
  {
    call <- match.call()
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    Cox <- (family$family == "Cox")
    if (Cox)
      require(survival)
    if (missing(data)) {
      data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data",  "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ",
                                                               method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm))
        names(Y) <- nm
    }
    # Interaction term, may require reordering
    if(length(grep(":", var))) {
      vars <- unlist(strsplit(var, split=":"))
      o <- order(match(vars, as.character(attr(mt,"variables"))))
      var <- paste(vars[o[1]],vars[o[2]],sep=":")
    }
    #
    if(Cox)
      if(!is.Surv(Y)){stop("for family 'Cox' the response variable has to be a survival object")}
    X <- if (!is.empty.model(mt))
      model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
      stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
      stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
      if (length(offset) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    if(!(paste(var) %in% colnames(X))) stop("var not a covariate in the formular")
    X <- cbind(X, rep(NA, nrow(X)))
    colnames(X)[ncol(X)] <- "resid"
    X[,"resid"] <- lm.fit(x = X[, -which(colnames(X) %in% c(paste(var),"resid")), drop = FALSE], y=X[, paste(var), drop = FALSE])$residuals
    ### Original data
    if(Cox){
      fit1 <- try(coxph.fit(x=X[, -which(colnames(X)==paste(var)), drop = FALSE], y=Y, strata=NULL, control=coxph.control(), method="efron", rownames=NULL) )
      fit2 <- try(coxph.fit(x=X[, -which(colnames(X) %in% c(paste(var),"resid")), drop = FALSE], y=Y, strata=NULL, control=coxph.control(), method="efron", rownames=NULL) )
      if(class(fit1)[1]=="try-error" | class(fit2)[1]=="try-error") return("error in fitting the Cox model for either the full or the reduced model")
      p.value.obs <- 1 - pchisq(abs(-2*fit1$loglik[2] + 2* fit2$loglik[2]), 1)
    } else {
      fit1 <- glm.fit(x = X[, -which(colnames(X)==paste(var)), drop = FALSE], y = Y, weights = weights, start = start,
                      etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt,
                                                                                                                                    "intercept") > 0)
      fit2 <- glm.fit(x = X[, -which(colnames(X) %in% c(paste(var),"resid")), drop = FALSE], y = Y, weights = weights, start = start,
                      etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt,
                                                                                                                                    "intercept") > 0)
      ###  Dispersion factor for model fit1
      df.r <- fit1$df.residual
      if(df.r > 0){
        if(any(fit1$weights == 0)){warning("observations with zero weight not used for calculating dispersion")}
        dispersion <- sum((fit1$weights * fit1$residuals^2)[fit1$weights > 0])/df.r
      }
      if(df.r==0){ dispersion <- NaN}
      if(df.r==0){warning("dispersion is Na")}
      if(fit1$family$family=="binomial" & dispersion>1.5){warning("estimated dispersion is > 1.5, rather use family = quasibinomial")}
      if(fit1$family$family=="binomial" & dispersion<0.5){warning("estimated dispersion is < 0.5, rather use family = quasibinomial")}
      if(fit1$family$family=="poisson" & dispersion>1.5){warning("estimated dispersion is > 1.5, rather use family = quasipoisson")}
      if(fit1$family$family=="poisson" & dispersion<0.5){warning("estimated dispersion is < 0.5, rather use family = quasipoisson")}
      estimated.dispersion <- dispersion
      if(fit1$family$family %in% c("poisson", "binomial")){dispersion <- 1}
      ###
      p.value.obs <- 1 - pchisq(abs(fit1$deviance - fit2$deviance)/dispersion, 1)
    }
    ### Permutations
    set.seed(seed)
    #
    if(Cox){
      devi  <- rep(NA, times=nrep)
      options(warn = -1)
      oldtime <- proc.time()[1]
      for (i in 1:nrep){devi[i] <- try(glm.perm(Y, X[, -which(colnames(X)==paste(var)), drop = FALSE], Family=family))}
      if(!Silent){print(c("execution time in minutes", round((proc.time()[1] - oldtime)/60, 2)))}
      options(warn = 0)
      psim <- 1 - pchisq(abs(devi + 2* fit2$loglik[2]), 1)
    } else {
      devi.disp <- matrix(0, ncol=2, nrow=nrep)
      options(warn = -1)
      oldtime <- proc.time()[1]
      for (i in 1:nrep){devi.disp[i,] <- glm.perm(Y, X[, -which(colnames(X)==paste(var)), drop = FALSE], Family=family)}
      if(!Silent){print(c("execution time in minutes", round((proc.time()[1] - oldtime)/60, 2)))}
      options(warn = 0)
      psim <- 1 - pchisq(abs(devi.disp[,1] - fit2$deviance)/devi.disp[,2], 1)
    }
    ### output of prr.test by Potter
    nrep.true <- sum(!is.na(psim))
    ret.val <- list(nobs = nrow(X), p0 = length(psim[psim <= p.value.obs])/nrep.true,
                    p005 = length(psim[psim <= 1.005 * p.value.obs])/nrep.true,
                    p01 = length(psim[psim <= 1.01 * p.value.obs])/nrep.true,
                    p02 = length(psim[psim <= 1.02 * p.value.obs])/nrep.true,
                    p04 = length(psim[psim <= 1.04 * p.value.obs])/nrep.true)
    names(ret.val$nobs) <- "number of observations used"
    names(ret.val$p0) <- "permutation p-value for simulated p-values <= observed p-value"
    names(ret.val$p005) <- "permutation p-value for simulated p-values <= 1.005 observed p-value"
    names(ret.val$p01) <- "permutation p-value for simulated p-values <= 1.01 observed p-value"
    names(ret.val$p02) <- "permutation p-value for simulated p-values <= 1.02 observed p-value"
    names(ret.val$p04) <- "permutation p-value for simulated p-values <= 1.04 observed p-value"
    ### new standard error output
    ret.stderr <- list(se.p0 = sqrt(ret.val$p0*(1-ret.val$p0)/nrep.true),
                       se.p005 = sqrt(ret.val$p005*(1-ret.val$p005)/nrep.true),
                       se.p01 = sqrt(ret.val$p01*(1-ret.val$p01)/nrep.true),
                       se.p02 = sqrt(ret.val$p02*(1-ret.val$p02)/nrep.true),
                       se.p04 = sqrt(ret.val$p04*(1-ret.val$p04)/nrep.true))
    names(ret.stderr$se.p0) <- NULL
    names(ret.stderr$se.p005) <- NULL
    names(ret.stderr$se.p01) <- NULL
    names(ret.stderr$se.p02) <- NULL
    names(ret.stderr$se.p04) <- NULL
    ### new output
    if (model)
      fit1$model <- mf
    fit1$na.action <- attr(mf, "na.action")
    if (x)
      fit1$x <- X
    if (!y)
      fit1$y <- NULL
    ###
    if(Cox){
      out <- c(list(fit1 = c(fit1,list(terms = mt, offset = offset, control = control, method = method,
                                       contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,mf))), fit2 = fit2, call = call, formula = formula,
                    seed=seed, fit1deviance=-2* fit2$loglik[2], fit2deviance=-2* fit2$loglik[2], LRstat= abs(-2* fit1$loglik[2] + 2* fit2$loglik[2]),  p.value.obs=p.value.obs, p.value.perm = ret.val, p.value.perm.se= ret.stderr, nobs= ret.val$nobs, var=var, nrep.true=nrep.true))
    } else {
      out <- c(list(fit1 = c(fit1,list(terms = mt, offset = offset, control = control, method = method,
                                       contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,mf))), fit2 = fit2, call = call, formula = formula,
                    seed=seed, fit1deviance=fit1$deviance, fit2deviance=fit2$deviance, Dispersion=dispersion, estimated.Dispersion = estimated.dispersion,
                    LRstat= abs(fit1$deviance - fit2$deviance)/dispersion,  p.value.obs=p.value.obs, p.value.perm = ret.val, p.value.perm.se= ret.stderr, nobs= ret.val$nobs, var=var, nrep.true=nrep.true))
    }
    class(out) <- "prr.test"
    return(out)
  }

#
# This is the new implementation for an arbitrary log-likelihood model
# Source code uses some code from these packages (GPL-2)
# https://cran.r-project.org/web/packages/pscl/index.html
# https://cran.r-project.org/web/packages/ZIBBSeqDiscovery/index.html
#

glm.fit.fixed <- function (x, y, ...) {
  fit <- glm.fit(x = x, y = y, ...)
  class(fit) <- c("glm", "lm")
  fit$loglik <- logLik(fit)[1]
  return(fit)
}

prr.test.new <- function (formula, data, subset, weights, na.action, offset, contrasts=NULL, test.var=NULL, n.rep=100, method="glm.fit.fixed", ...)
{
  cal <- match.call()
  #if (is.character(family))
  #    family <- get(family, mode = "function", envir = parent.frame())
  #if (is.function(family))
  #    family <- family()
  #if (is.null(family$family)) {
  #    print(family)
  #    stop("'family' not recognized")
  #}
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  #if (identical(method, "glm.fit"))
  #    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  # Which columns in X are related to test.var?
  term.labels  <- attr(mt, "term.labels")
  term.columns <- attr(X, "assign")
  if (attr(mt, "intercept")) {
    term.labels  <- c("(Intercept)", term.labels)
    term.columns <- term.columns + 1
  }
  cols.testvar <- term.labels[term.columns] == test.var
  # Predict test.var from other variables
  fit <- lm.fit(X[,!cols.testvar], X[,cols.testvar])
  R   <- fit$residual
  if (is.null(dim(R)) | (length(dim(R)) == 1L))
    dim(R) <- c(length(R), 1)

  # Calculate observed p-value
  X0 <- X[,!cols.testvar]
  XR <- cbind(X0, R)
  f     <- call(if (is.function(method)) "method" else method)
  fit.1 <- eval(as.call(c(as.list(f), list(x = XR, y = Y, weights = weights, offset = offset, intercept = attr(mt, "intercept") > 0),
                          substitute(...()))))
  fit.2 <- eval(as.call(c(as.list(f), list(x = X0, y = Y, weights = weights, offset = offset, intercept = attr(mt, "intercept") > 0),
                          substitute(...()))))
  ll.1  <- fit.1$loglik
  ll.2  <- fit.2$loglik
  lr.df <- fit.2$df.residual - fit.1$df.residual
  p.value.obs <- pchisq(2 * (ll.1 - ll.2), lr.df, lower.tail=F)

  # Calculate simulated p-value
  ll.r <- rep(0, length(n.rep))
  for (i in 1:n.rep) {
    rnd <- sample(1:nrow(X), replace = FALSE)
    XR  <- cbind(X0, R[rnd,])
    fit.r <- eval(as.call(c(as.list(f), list(x = XR, y = Y, weights = weights, offset = offset, intercept = attr(mt, "intercept") > 0),
                            substitute(...()))))
    ll.r[i] <- fit.r$loglik
  }
  p.value.sim <- mean(ll.1 < ll.r)
  return(list(p.value.obs=p.value.obs, p.value.sim=p.value.sim))
}


model_offset_2 <- function(x, terms = NULL, offset = TRUE)
{
  if(is.null(terms)) terms <- attr(x, "terms")
  offsets <- attr(terms, "offset")
  if(length(offsets) > 0) {
    ans <- if(offset) x$"(offset)" else NULL
    if(is.null(ans)) ans <- 0
    for(i in offsets) ans <- ans + x[[deparse(attr(terms, "variables")[[i + 1]])]]
    ans
  }
  else {
    ans <- if(offset) x$"(offset)" else NULL
  }
  if(!is.null(ans) && !is.numeric(ans)) stop("'offset' must be numeric")
  ans
}

linkobj = make.link("logit") #linkinv= function(x) 1 / (1 + exp(-x))

#Negative Binomial likelihood P(N=n)
NegBin <- function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
  kx <- ncol(X)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  theta <- exp(parms[kx + 1])
  loglik <- sum(weights * suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE)))
  return(loglik)
}

# Beta-Binomial likelihood P(N=n)
BetaBin <- function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  theta <- linkobj$linkinv(parms[kx + 1])
  s1 <- as.vector(mu * (1 - theta)/theta)
  s2 <- as.vector((1 - theta)/theta * (1 - mu))
  loglike <- sum(weights * (lchoose(y.c, y) +  lbeta(y + s1, y.c - y + s2) - lbeta(s1, s2)))
  return(loglike)
}

#Poisson zero-inflated likelihood P(N=n)
ziPoisson <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y0 <- Y <= 0
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkobj$linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  loglik0 <- log(phi + exp(log(1 - phi) - mu))
  loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
  loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
  return(loglik)
}
#Gradient
gradPoisson <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  eta <- as.vector(X %*% parms[1:kx] + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
  muz <- linkobj$linkinv(etaz)
  clogdens0 <- -mu
  dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)
  wres_count <- ifelse(Y1, Y - mu, -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(mu)))
  wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
                      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
  grad <- colSums(cbind(wres_count * weights * X, wres_zero * weights * Z))
  return(grad)
}

#Negative Binomial zero-inflated likelihood P(N=n)
ziNegBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y0 <- Y <= 0
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkobj$linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  theta <- exp(parms[(kx + kz) + 1])
  loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
  return(loglik)
}
# Gradient
gradNegBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  eta <- as.vector(X %*% parms[1:kx] + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
  muz <- linkobj$linkinv(etaz)
  theta <- exp(parms[(kx + kz) + 1])
  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +  clogdens0)
  wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta),
                       -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) -
                              log(mu + theta) + log(mu)))
  wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
                      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
  wres_theta <- theta * ifelse(Y1, digamma(Y + theta) -
                                 digamma(theta) + log(theta) - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 -
                                                                                   muz) + clogdens0) * (log(theta) - log(mu + theta) +
                                                                                                          1 - theta/(mu + theta)))
  grad <- colSums(cbind(wres_count * weights * X, wres_zero * weights *
                          Z, wres_theta))
  return(grad)
}

# Binomial zero-inflated likelihood P(N=n)
ziBinomial <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  kz <- NCOL(Z)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  phi <- linkobj$linkinv(as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  Y0 <- y == 0
  loglike0 <- log(phi + exp(log(1 - phi) + y.c*log(1 - mu)))
  loglike1 <- log(1 - phi) + lchoose(y.c, y) +  y*log(mu) + (y.c - y)*log(1 - mu)
  loglike <- sum(weights[Y0] * loglike0[Y0]) + sum(weights[!Y0] * loglike1[!Y0])
  return(loglike)
}

# Beta-Binomial zero-inflated likelihood P(N=n)
ziBetaBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  kz <- NCOL(Z)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  phi <- linkobj$linkinv(as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  theta <- linkobj$linkinv(parms[kx + kz + 1])
  s1 <- as.vector(mu * (1 - theta)/theta)
  s2 <- as.vector((1 - theta)/theta * (1 - mu))
  Y0 <- y == 0
  loglike0 <- log(phi + exp(log(1 - phi) + lbeta(s1, y.c + s2) - lbeta(s1, s2)))
  loglike1 <- log(1 - phi) + lchoose(y.c, y) +  lbeta(y + s1, y.c - y + s2) - lbeta(s1, s2)
  loglike <- sum(weights[Y0] * loglike0[Y0]) + sum(weights[!Y0] * loglike1[!Y0])
  return(loglike)
}

fitdist.control <- function (method = "BFGS", maxit = 10000, trace = FALSE, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, trace = trace, start = start)
  rval <- c(rval, list(...))
  rval$fnscale <- -1
  rval$hessian <- TRUE
  if (is.null(rval$reltol))
    rval$reltol <- .Machine$double.eps^(1/1.6)
  return(rval)
}

fitdist <- function (X, Y, Z=NULL, offsetx=NULL, offsetz=NULL, weights=NULL,
                     dist = c("negbin", "betabin", "zipoisson", "zinegbin", "zibinomial", "zibetabin"),
                     control = fitdist.control(...), ...)
{
  dist <- match.arg(dist)
  #TODO: implement gradient for binomial / beta-binomial
  loglikfun <- switch(dist, negbin=NegBin, betabin=BetaBin, zipoisson = ziPoisson, zinegbin = ziNegBin, zibinomial=ziBinomial, zibetabin = ziBetaBin)
  gradfun <- switch(dist, negbin=NULL, betabin=NULL, zipoisson = gradPoisson, zinegbin = gradNegBin, zibinomial=NULL, zibetabin = NULL)

  zero.inflated <- startsWith(dist, "zi")
  over.disperse <- dist %in% c("negbin", "betabin", "zinegbin", "zibetabin")

  if (control$trace) {
    if (zero.inflated)
      cat("Zero-inflated Count Model\n", dist, "\n", sep = "")
    else
      cat("Count Model\n", dist, "\n", sep = "")
  }

  n <- NROW(Y)
  kx <- NCOL(X)
  kz <- ifelse(zero.inflated, NCOL(Z), 0)

  if (is.null(weights))
    weights = rep(1, n)
  if (is.null(offsetx))
    offsetx = rep(0,n)
  if (zero.inflated & is.null(offsetz))
    offsetz = rep(0,n)

  start <- control$start
  if (is.null(start)) {
    if (control$trace)
      cat("generating starting values...")
    if (dist %in% c("negbin", "zipoisson", "zinegbin")) {
      if (zero.inflated)
        model_zero <- glm.fit(Z, as.integer(Y == 0), weights = weights, family = binomial(link = "logit"), offset = offsetz)
      else
        model_zero <- list(coefficients=NULL)
      model_count <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)
      theta <- 1
    } else if (dist %in% c("betabin", "zibinomial", "zibetabin")) {
      y <- as.vector(Y[,1])
      y.c <- as.vector(rowSums(Y))
      phat <- y/y.c
      phat[phat <= 0] <- 1/(2 * y.c[phat == 0])
      phat[phat >= 1] <- 1 - 1/(2 * y.c[phat == 1])
      if (zero.inflated)
        model_zero <- lm.fit(Z, as.numeric(y == 0))
      else
        model_zero <- list(coefficients=NULL)
      model_count <- lm.fit(X, log(phat/(1 - phat)))
      model.p <- mean(y) / mean(y.c)
      phihat <- (1/(mean(y.c) - 1)) * (var(y)/(mean(y.c) * model.p * (1 - model.p)) - 1)
      theta <- phihat/(1 - phihat)
    }
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients,
                  theta = if (over.disperse) log(theta) else NULL)
  }

  method <- control$method
  hessian <- control$hessian
  control$method <- control$hessian <- control$start <- NULL
  if (control$trace)
    cat("calling optim() for ML estimation:\n")
  fit <- optim(fn = loglikfun, gr = gradfun, par = c(start$count, start$zero, start$theta),
               X = X, Z = Z, Y = Y, offsetx = offsetx, offsetz = offsetz, weights = weights,
               method = method, hessian = hessian, control = control)
  if (fit$convergence > 0)
    warning("optimization failed to converge")
  coefc <- fit$par[1:kx]
  names(coefc) <- names(start$count) <- colnames(X)
  if (zero.inflated) {
    coefz <- fit$par[(kx + 1):(kx + kz)]
    names(coefz) <- names(start$zero) <- colnames(Z)
  } else {
    coefz <- NULL
  }
  if (over.disperse)
    theta <- fit$par[(kx + kz + 1)]
  else
    theta <- NULL
  vc <- tryCatch(-solve(as.matrix(fit$hessian)), error = function(e) {
    warning(e$message, call = FALSE)
    k <- nrow(as.matrix(fit$hessian))
    return(matrix(NA, k, k))
  })
  cols.count <- paste("count", colnames(X), sep = "_")
  cols.zero <- if (zero.inflated) paste("zero", colnames(Z), sep = "_") else NULL
  cols.theta <- if (over.disperse) "theta" else NULL
  colnames(vc) <- rownames(vc) <- c(cols.count, cols.zero, cols.theta)
  mu <- exp(X %*% coefc + offsetx)[, 1]
  if (zero.inflated) {
    phi <- linkobj$linkinv(Z %*% coefz + offsetz)[, 1]
    Yhat <- (1 - phi) * mu
  } else {
    Yhat <- mu
  }
  res <- sqrt(weights) * (Y - Yhat)
  nobs <- sum(weights > 0)
  rval <- list(
    coefficients = list(count = coefc, zero = coefz, theta=theta), vcov = vc,
    residuals = res, fitted.values = Yhat, dist = dist, loglik = fit$value, n = nobs,
    df.null = nobs - 2, df.residual = nobs - (kx + kz + (dist %in% c("negbin", "betabin"))),
    converged = fit$convergence < 1,
    optim = fit, start = start,
    weights = if (identical(weights, rep.int(1, n))) NULL else weights,
    offset = list(count = if (identical(offsetx, rep.int(0, n))) NULL else offsetx,
                  zero = if (identical(offsetz, rep.int(0, n))) NULL else offsetz),
    contrasts = list(count = attr(X, "contrasts"),
                     zero = if(zero.inflated) attr(Z, "contrasts") else NULL))
  #class(rval) <- "zeroinfl"
  return(rval)
}

prr.test.new.zeroinfl <- function (formula, data, subset, weights, na.action, offset,
                                   test.var = NULL, test = "both", n.rep=100, dist="zipoisson",
                                   trace=F, xyz=F, ...)
{
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset",
               "na.action", "weights", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  # Extract formula y ~ a + b + c | d + e into ffc = y ~ a + b + c, ffz = ~ d + e
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
    ff <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
    formula[[3]][1] <- call("+")
    mf$formula <- formula
  }
  else {
    ffc <- ff <- formula
    ffz <- NULL
    #ffz[[2]] <- NULL
  }

  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  # Response
  Y <- model.response(mf, "numeric")
  # Weights
  n <- NROW(Y)
  weights <- model.weights(mf)
  if (is.null(weights))
    weights <- 1
  if (length(weights) == 1)
    weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  # X: count model matrix
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
  if (is.null(offsetx))
    offsetx <- 0
  if (length(offsetx) == 1)
    offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  # Which columns in X are related to test.var?
  X.term.labels  <- attr(mtX, "term.labels")
  X.term.columns <- attr(X, "assign")
  if (attr(mtX, "intercept")) {
    X.term.labels  <- c("(Intercept)", X.term.labels)
    X.term.columns <- X.term.columns + 1
  }

  # Z: zero inflation model matrix
  if(!is.null(ffz)) {
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~.), data = data)
    Z <- model.matrix(mtZ, mf)
    offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
    if (is.null(offsetz))
      offsetz <- 0
    if (length(offsetz) == 1)
      offsetz <- rep.int(offsetz, n)
    offsetz <- as.vector(offsetz)
    # Which columns in X are related to test.var?
    Z.term.labels  <- attr(mtZ, "term.labels")
    Z.term.columns <- attr(Z, "assign")
    if (attr(mtZ, "intercept")) {
      Z.term.labels  <- c("(Intercept)", Z.term.labels)
      Z.term.columns <- Z.term.columns + 1
    }
  } else {
    Z <- NULL
    offsetz <- NULL
  }

  X.cols.testvar = X.term.labels[X.term.columns] == test.var
  Z.cols.testvar = if (is.null(ffz)) NULL else Z.term.labels[Z.term.columns] == test.var

  # Check if we should only test the "count" or "zero" component
  if (test == "count")
    Z.cols.testvar = NULL
  if (test == "zero")
    X.cols.testvar = NULL

  # Required count component
  if (!is.null(X.cols.testvar)) {

    # Predict test.var from other variables
    fit.X <- lm.fit(X[,!X.cols.testvar], X[,X.cols.testvar])
    R.X   <- fit.X$residual
    if (is.null(dim(R.X)) | (length(dim(R.X)) == 1L))
      dim(R.X) <- c(length(R.X), 1)

    # Model matrix without variable of interest
    # and with residuals of variable of interest
    X0 <- X[,!X.cols.testvar]
    XR <- cbind(X0, R.X)
  } else {
    XR <- X0 <- X
  }

  # Possible zero-inflation component
  if (!is.null(Z.cols.testvar)) {
    # Predict test.var from other variables
    fit.Z <- lm.fit(Z[,!Z.cols.testvar], Z[,Z.cols.testvar])
    R.Z   <- fit.Z$residual
    if (is.null(dim(R.Z)) | (length(dim(R.Z)) == 1L))
      dim(R.Z) <- c(length(R.Z), 1)

    # Model matrix without variable of interest
    # and with residuals of variable of interest
    Z0 <- Z[,!Z.cols.testvar]
    ZR <- cbind(Z0, R.Z)
  } else {
    ZR <- Z0 <- Z
  }

  # Calculate observed p-value
  fit.1 <- fitdist(XR, Y, ZR, offsetx, offsetz, weights, dist=dist, control=fitdist.control(trace=trace))
  fit.2 <- fitdist(X0, Y, Z0, offsetx, offsetz, weights, dist=dist, control=fitdist.control(trace=trace))
  ll.1    <- fit.1$loglik
  ll.2    <- fit.2$loglik
  lr.df <- fit.2$df.residual - fit.1$df.residual
  p.value.obs <- pchisq(2 * (ll.1 - ll.2), lr.df, lower.tail=F)

  # Calculate simulated p-value
  ll.r   <- rep(0, length(n.rep))
  for (i in 1:n.rep) {
    rnd <- sample(1:nrow(X), replace = FALSE)
    if (!is.null(X.cols.testvar))
      XR <- cbind(X0, R.X[rnd,])
    if (!is.null(Z.cols.testvar))
      ZR <- cbind(Z0, R.Z[rnd,])
    fit.r   <- fitdist(XR, Y, ZR, offsetx, offsetz, weights, dist=dist, control=fitdist.control(trace=trace))
    ll.r[i]   <- fit.r$loglik
  }
  p.value.sim   <- mean(ll.1 < ll.r)

  results <- list(p.value.obs=p.value.obs, p.value.sim=p.value.sim)

  if (xyz)
    results$data <- list(X=X, Z=Z, Y=Y, offsetx=offsetx, offsetz=offsetz, weights=weights)

  return(results)
}
