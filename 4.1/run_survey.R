.libPaths("~/R/hpc_libs")

library(survey)
library(stats)

# from postlink v0.1.2
constrained_logistic_regression <- function(X, y, bound, cmax.iter) {
 X <- as.matrix(X)
 d <- ncol(X)
 barx <- colMeans(X)
 tau <- 0.1
 
 nloglik_logistic <- function(mu) {
  cutoff <- 1E-10
  flag0 <- (mu < cutoff)
  flag1 <- (mu > 1 - cutoff)
  flag <- flag0 | flag1
  nloglik_nonflag <- sum(-(y[!flag] * log(mu[!flag]) + (1 - y[!flag]) * log(1 - mu[!flag])))
  
  if (any(flag0) && max(y[flag0]) > cutoff) {
   warning("Fitted probabilities zero with non-matching response; evaluation may be inaccurate.")
  }
  if (any(flag1) && min(y[flag1]) < 1 - cutoff) {
   warning("Fitted probabilities one with non-matching response; evaluation may be inaccurate.")
  }
  return(nloglik_nonflag)
 }
 
 grad_Hess <- function(mu) {
  grad <- -crossprod(X, y - mu)
  w <- sqrt(mu * (1 - mu))
  Xw <- sweep(X, MARGIN = 1, STATS = w, FUN = "*")
  Hess <- crossprod(Xw)
  return(list(grad = grad, Hess = Hess))
 }
 
 betacur <- c(bound, rep(0, d - 1))
 mu <- plogis(X %*% betacur)
 iter <- 1
 maxiter <- cmax.iter
 tol <- 1E-8
 objs <- numeric(maxiter)
 objs[1] <- nloglik_logistic(mu)
 
 while (iter < maxiter) {
  gH <- grad_Hess(mu)
  gr <- gH$grad
  H <- gH$Hess
  rhs <- -gr + (H %*% betacur)
  betanew <- solve(H, rhs)
  eta_tmp <- X %*% betanew
  
  if (mean(eta_tmp) > bound) {
   betanew_theta <- solve(rbind(cbind(H, barx), c(barx, 0)), c(rhs, bound))
   betanew <- betanew_theta[1:d]
  }
  
  m <- 0
  upd <- betanew - betacur
  munew <- plogis(X %*% betanew)
  objnew <- nloglik_logistic(munew)
  
  while ((objs[iter] - objnew) < tau * (1/2)^m * sum(upd * (-gr))) {
   m <- m + 1
   if ((1/2)^m < tol^2) {
    betanew <- betacur
    break
   } else {
    betanew <- betacur + upd * (1/2)^m
    munew <- plogis(X %*% betanew)
    objnew <- nloglik_logistic(munew)
   }
  }
  betacur <- betanew
  mu <- munew
  iter <- iter + 1
  objs[iter] <- objnew
  if (iter >= 2 && objs[iter - 1] - objs[iter] < tol) break
 }
 return(list(beta = betacur, objs = objs))
}

glmMixture <- function(x, y, family,
                       z, m.rate = NULL, safe.matches = NULL,
                       control = list(), ...) {
 
 # 1. Set-up and Arguments
 dots <- list(...)
 con <- list(init.beta = NULL, init.gamma = NULL, fy = NULL,
             max.iter = 1000, tol = 1E-4, cmax.iter = 1000)
 con[names(control)] <- control
 con[names(dots)] <- dots
 
 x <- as.matrix(x)
 y <- as.numeric(y)
 n <- nrow(x)
 p <- ncol(x)
 d <- p
 z <- as.matrix(z)
 
 if (is.null(safe.matches)) safe.matches <- rep(FALSE, n)
 
 if (!is.null(m.rate)) {
  logitbound <- -log((1 - m.rate) / m.rate)
 } else {
  logitbound <- NULL
 }
 
 family_arg <- if (!is.null(family)) family else stats::gaussian()
 if (is.character(family_arg)) {
  if (tolower(family_arg) == "gamma") family <- stats::Gamma(link = "log")
  else {
   family_fun <- try(get(family_arg, mode = "function", envir = parent.frame()), silent = TRUE)
   if (inherits(family_fun, "try-error")) stop(paste("Invalid family argument:", family_arg))
   family <- family_fun()
  }
 } else if (is.function(family_arg)) {
  family <- family_arg()
 } else {
  family <- family_arg
 }
 family_name <- tolower(family$family)
 
 if (family_name == "gaussian" && family$link != "identity") {
  stop("The 'gaussian' family is currently restricted to the 'identity' link.")
 }
 
 # Marginal Density (fy) Estimation
 fy <- con$fy
 if (is.null(fy)) {
  if (family_name == "gaussian") {
   fy <- stats::dnorm(y, mean = mean(y), sd = stats::sd(y))
  } else if (family_name == "binomial") {
   fy <- (mean(y)^y) * (1 - mean(y))^(1 - y)
  } else {
   kde_obj <- stats::density(y)
   g <- stats::approxfun(x = kde_obj$x, y = kde_obj$y, method = "linear")
   fy <- g(y)
  }
 }
 # prevent 0s in custom marginals or KDE boundaries
 if (!is.null(fy)) {
  fy[is.na(fy) | fy <= 0] <- 1e-10
 }
 
 # Initialization: Beta & Shape
 beta_cur <- con$init.beta
 
 # Define shape_update globally so it is always available for the EM loop
 shape_update <- function(mu_val, mme, pcur) {
  f <- function(shape) {
   -sum(pcur * stats::dgamma(y, shape = shape, rate = shape / mu_val, log = TRUE))
  }
  lower_bound <- max(1e-4, 0.1 * mme)
  upper_bound <- max(10, 10 * mme)
  stats::optimize(f = f, lower = lower_bound, upper = upper_bound)$minimum
 }
 
 if (is.null(beta_cur)) {
  if (family_name != "gamma") {
   init_fit <- stats::glm.fit(x, y, family = family)
   beta_cur <- stats::coef(init_fit)
  } else {
   init_fit <- stats::glm.fit(x, y, family = stats::Gamma(link = "log"))
   beta_cur <- stats::coef(init_fit)
  }
 }
 
 # Always initialize dispersion parameters, regardless of whether beta was provided
 shape_cur <- NA
 stdcur <- NA
 
 if (family_name == "gaussian") {
  stdcur <- sqrt(sum((y - x %*% beta_cur)^2) / max(1, n - p))
 } else if (family_name == "gamma") {
  mu_init <- pmax(family$linkinv(x %*% beta_cur), 1e-10)
  disp_init <- sum(((y - mu_init)/mu_init)^2) / max(1, n - p)
  if (is.na(disp_init) || disp_init <= 0) disp_init <- 1
  shape_cur <- shape_update(mu_init, mme = 1/disp_init, rep(1, n))
 }
 
 # Initialization: Gamma
 gamma_cur <- con$init.gamma
 if (!is.null(gamma_cur)) {
  gamma_cur <- -1 * gamma_cur
 }
 if (is.null(gamma_cur)) {
  if (!is.null(logitbound)) {
   if (ncol(z) == 1 && all(z == 1)) {
    gamma_cur <- rep(-logitbound, ncol(z))
   } else {
    gamma_cur <- c(max(-logitbound, 0), rep(0, ncol(z) - 1))
   }
  } else {
   gamma_cur <- rep(0, ncol(z))
  }
 }
 
 # Helpers
 hgamma <- function(eta) {
  pi <- stats::plogis(eta)
  list(fun = pi, dfun = pi * (1 - pi), d2fun = pi * (1 - pi) * (1 - 2 * pi))
 }
 
 fymu_eval <- function(mu, sigma, sub, shape) {
  mu_safe <- mu[sub]
  if (family_name %in% c("poisson", "gamma")) {
   mu_safe <- pmax(mu_safe, 1e-10)
  } else if (family_name == "binomial") {
   mu_safe <- pmax(pmin(mu_safe, 1 - 1e-10), 1e-10)
  }
  
  if (family_name == "gaussian") return(stats::dnorm((y[sub] - mu_safe), sd = sigma))
  if (family_name == "poisson") return(stats::dpois(y[sub], mu_safe))
  if (family_name == "binomial") return(stats::dbinom(y[sub], 1, mu_safe))
  if (family_name == "gamma") return(stats::dgamma(y[sub], shape, shape / mu_safe))
 }
 
 nloglik <- function(mu, sigma, hs, shape) {
  term1 <- hs[!safe.matches] * fymu_eval(mu, sigma, !safe.matches, shape) + (1 - hs[!safe.matches]) * fy[!safe.matches]
  term2 <- fymu_eval(mu, sigma, safe.matches, shape)
  
  sum(-log(pmax(term1, 1e-300))) - sum(log(pmax(term2, 1e-300)))
 }
 
 # 2. EM-Algorithm
 eta_cur <- x %*% beta_cur
 mu_cur <- family$linkinv(eta_cur)
 hs <- hgamma(z %*% gamma_cur)$fun
 hs[safe.matches] <- 1
 p_cur <- rep(0, n)
 
 iter <- 1
 objs <- numeric(con$max.iter)
 objs[iter] <- nloglik(mu_cur, stdcur, hs, shape_cur)
 
 converged_flag <- FALSE
 w_glm_fit <- NULL
 
 while (iter < con$max.iter) {
  # E-Step
  num <- hs[!safe.matches] * fymu_eval(mu_cur, stdcur, !safe.matches, shape_cur)
  denom <- num + (1 - hs[!safe.matches]) * fy[!safe.matches]
  denom[denom <= 0] <- 1e-10
  
  p_cur[!safe.matches] <- pmin(pmax(num / denom, 1e-10), 1 - 1e-10)
  p_cur[safe.matches] <- 1
  
  # M-Step: Gamma
  if (!is.null(logitbound)) {
   if (ncol(z) == 1 && all(z == 1)) {
    temp_glm_h <- tryCatch(stats::glm.fit(z[!safe.matches, , drop = FALSE], p_cur[!safe.matches], family = stats::quasibinomial()), error = function(e) NULL)
    if (!is.null(temp_glm_h) && !anyNA(stats::coef(temp_glm_h))) {
     gamma_cur <- max(stats::coef(temp_glm_h), -logitbound)
    } else {
     warning("Mismatch model failed to find valid coefficients. Terminating EM early.")
     break
    }
   } else {
    temp_glm_h <- constrained_logistic_regression(z[!safe.matches, , drop = FALSE], 1 - p_cur[!safe.matches], logitbound, con$cmax.iter)
    if (!is.null(temp_glm_h) && !anyNA(temp_glm_h$beta)) {
     gamma_cur <- -temp_glm_h$beta
    } else {
     warning("Mismatch model failed to find valid coefficients. Terminating EM early.")
     break
    }
   }
  } else {
   temp_glm_h <- tryCatch(stats::glm.fit(z[!safe.matches, , drop = FALSE], p_cur[!safe.matches], family = stats::quasibinomial()), error = function(e) NULL)
   if (!is.null(temp_glm_h) && !anyNA(stats::coef(temp_glm_h))) {
    gamma_cur <- stats::coef(temp_glm_h)
   } else {
    warning("Mismatch model failed to find valid coefficients. Terminating EM early.")
    break
   }
  }
  hs[!safe.matches] <- hgamma(z[!safe.matches, , drop = FALSE] %*% gamma_cur)$fun
  
  # M-Step: Beta & Dispersion
  if (family_name != "gamma") {
   temp_glm_fit <- tryCatch(stats::glm.fit(x, y, weights = p_cur, family = family, start = beta_cur), error = function(e) NULL)
   if (is.null(temp_glm_fit) || anyNA(stats::coef(temp_glm_fit))) {
    temp_glm_fit <- tryCatch(stats::glm.fit(x, y, weights = p_cur, family = family), error = function(e) NULL)
   }
   
   if (!is.null(temp_glm_fit) && !anyNA(stats::coef(temp_glm_fit))) {
    w_glm_fit <- temp_glm_fit
    beta_cur <- stats::coef(w_glm_fit)
   } else {
    warning("Outcome model failed to converge in M-step. Terminating EM algorithm early.")
    break
   }
  } else {
   temp_glm_fit <- tryCatch(stats::glm.fit(x, y, weights = p_cur, family = family, start = beta_cur), error = function(e) NULL)
   if (is.null(temp_glm_fit) || anyNA(stats::coef(temp_glm_fit))) {
    temp_glm_fit <- tryCatch(stats::glm.fit(x, y, weights = p_cur, family = family), error = function(e) NULL)
   }
   
   if (!is.null(temp_glm_fit) && !anyNA(stats::coef(temp_glm_fit))) {
    w_glm_fit <- temp_glm_fit
    beta_cur <- stats::coef(w_glm_fit)
    mu_w <- pmax(w_glm_fit$fitted.values, 1e-10)
    disp_w <- sum(p_cur * ((y - mu_w)/mu_w)^2) / max(sum(p_cur), 1)
    if (is.na(disp_w) || disp_w <= 0) disp_w <- 1
    
    try_shape <- tryCatch(shape_update(mu_w, mme = 1/disp_w, p_cur), error = function(e) NA)
    if (!is.na(try_shape)) shape_cur <- try_shape
   } else {
    warning("Outcome model failed to converge in M-step. Terminating EM algorithm early.")
    break
   }
  }
  
  eta_cur <- x %*% beta_cur
  mu_cur <- family$linkinv(eta_cur)
  
  if (family_name == "gaussian") {
   stdcur <- sqrt(stats::weighted.mean((y - mu_cur)^2, w = p_cur))
  }
  
  iter <- iter + 1
  objs[iter] <- nloglik(mu_cur, stdcur, hs, shape_cur)
  
  if (is.na(objs[iter])) {
   warning("EM algorithm did not converge. NA objective value occurred.")
   break
  }
  if (abs(objs[iter] - objs[iter - 1]) < con$tol) {
   converged_flag <- TRUE
   break
  }
 }
 
 # 3. Standard Errors
 hgamma_eval <- hgamma(z %*% as.matrix(gamma_cur))
 
 if (any(safe.matches)) {
  hgamma_eval$fun[safe.matches] <- 1
  hgamma_eval$dfun[safe.matches] <- 0
  hgamma_eval$d2fun[safe.matches] <- 0
 }
 
 if (family_name == "gaussian") {
  fymu_all_eval <- (function(mu, std) {
   fun <- stats::dnorm((y - mu), sd = std)
   d_fun_beta <- fun * (y - mu) / (std^2)
   d_fun_sigma <- fun * ((y - mu)^2 / std^3 - 1 / std)
   d2_fun_beta <- d_fun_beta * (y - mu) / (std^2) - fun / (std^2)
   d2_fun_sigma <- d_fun_sigma * ((y - mu)^2 / std^3 - 1 / std) + fun * (1 / std^2 - 3 * (y - mu)^2 / std^4)
   d2_fun_beta_sigma <- d_fun_sigma * (y - mu) / (std^2) - 2 * fun * (y - mu) / (std^3)
   list(fun = fun, dfun_beta = d_fun_beta, dfun_sigma = d_fun_sigma,
        d2fun_beta = d2_fun_beta, d2fun_sigma = d2_fun_sigma, d2fun_beta_sigma = d2_fun_beta_sigma)
  })(mu_cur, stdcur)
  
  mixprob <- pmax(fy * (1 - hgamma_eval$fun) + hgamma_eval$fun * fymu_all_eval$fun, 1e-10)
  
  w_beta_score <- (-1) * fymu_all_eval$dfun_beta * hgamma_eval$fun / mixprob
  w_sigma_score <- (-1) * fymu_all_eval$dfun_sigma * hgamma_eval$fun / mixprob
  w_gamma_score <- (-1) * (fymu_all_eval$fun - fy) * hgamma_eval$dfun / mixprob
  
  Xw1 <- sweep(x, 1, w_beta_score, "*")
  Xw2 <- matrix(w_sigma_score, n, 1)
  Deltaw3 <- sweep(z, 1, w_gamma_score, "*")
  meat <- crossprod(cbind(Xw1, Xw2, Deltaw3))
  
  w_beta2_hess <- (-(hgamma_eval$fun * fymu_all_eval$d2fun_beta) / mixprob) + w_beta_score^2
  w_sigma2_hess <- (-(hgamma_eval$fun * fymu_all_eval$d2fun_sigma) / mixprob) + w_sigma_score^2
  w_gamma2_hess <- ((-(fymu_all_eval$fun - fy) * hgamma_eval$d2fun) / mixprob) + w_gamma_score^2
  
  w_beta_gamma_hess <- (-(fymu_all_eval$dfun_beta * hgamma_eval$dfun) / mixprob) +
   ((fymu_all_eval$fun - fy) * hgamma_eval$fun * hgamma_eval$dfun * fymu_all_eval$dfun_beta) / (mixprob^2)
  w_sigma_gamma_hess <- (-(fymu_all_eval$dfun_sigma * hgamma_eval$dfun) / mixprob) +
   ((fymu_all_eval$fun - fy) * hgamma_eval$fun * hgamma_eval$dfun * fymu_all_eval$dfun_sigma) / (mixprob^2)
  
  w_beta_sigma_hess <- (-(fymu_all_eval$d2fun_beta_sigma * hgamma_eval$fun) / mixprob) + w_beta_score * w_sigma_score
  
  Xw4 <- sweep(x, 1, w_beta2_hess, "*")
  Deltaw6 <- sweep(z, 1, w_gamma2_hess, "*")
  Xw5 <- sweep(x, 1, w_beta_gamma_hess, "*")
  
  Hess <- matrix(0, d + 1 + ncol(z), d + 1 + ncol(z))
  one_v <- matrix(1, n, 1)
  
  Hess[1:d, 1:d] <- crossprod(x, Xw4)
  Hess[d+1, d+1] <- crossprod(one_v, sweep(one_v, 1, w_sigma2_hess, "*"))
  Hess[(d+2):ncol(Hess), (d+2):ncol(Hess)] <- crossprod(z, Deltaw6)
  
  Hess[1:(d+1), (d+2):ncol(Hess)] <- rbind(crossprod(Xw5, z),
                                           crossprod(sweep(one_v, 1, w_sigma_gamma_hess, "*"), z))
  Hess[(d+2):ncol(Hess), 1:(d+1)] <- t(Hess[1:(d+1), (d+2):ncol(Hess)])
  
  Hess[1:d, d+1] <- crossprod(x, sweep(one_v, 1, w_beta_sigma_hess, "*"))
  Hess[d+1, 1:d] <- t(Hess[1:d, d+1])
  
 } else {
  fymu_all_eval <- (function(eta, family_obj, shape) {
   mu_sub <- family_obj$linkinv(eta)
   
   if (family_obj$family %in% c("poisson", "Gamma")) mu_sub <- pmax(mu_sub, 1e-10)
   if (family_obj$family == "binomial") mu_sub <- pmax(pmin(mu_sub, 1 - 1e-10), 1e-10)
   
   if (family_obj$family == "poisson") {
    fun <- stats::dpois(y, mu_sub)
   } else if (family_obj$family == "binomial") {
    fun <- stats::dbinom(y, 1, mu_sub)
   } else if (family_obj$family == "Gamma") {
    fun <- stats::dgamma(y, shape, shape / mu_sub)
   } else {
    fun <- rep(1, length(y))
   }
   
   mu_prime <- family_obj$mu.eta(eta)
   var_mu <- family_obj$variance(mu_sub)
   
   if (family_obj$family == "Gamma") {
    score_eta <- (y - mu_sub) / var_mu * mu_prime * shape
    expected_hess <- - (mu_prime^2) / var_mu * shape
   } else {
    score_eta <- (y - mu_sub) / var_mu * mu_prime
    expected_hess <- - (mu_prime^2) / var_mu
   }
   
   dfun <- fun * score_eta
   d2fun <- fun * (score_eta^2 + expected_hess)
   
   list(fun = fun, dfun = dfun, d2fun = d2fun)
  })(eta_cur, family, shape_cur)
  
  mixprob <- pmax(fy * (1 - hgamma_eval$fun) + hgamma_eval$fun * fymu_all_eval$fun, 1e-10)
  
  w_beta_score <- (-1) * fymu_all_eval$dfun * hgamma_eval$fun / mixprob
  w_gamma_score <- (-1) * (fymu_all_eval$fun - fy) * hgamma_eval$dfun / mixprob
  
  Xw1 <- sweep(x, 1, w_beta_score, "*")
  Deltaw3 <- sweep(z, 1, w_gamma_score, "*")
  meat <- crossprod(cbind(Xw1, Deltaw3))
  
  w_beta2_hess <- (-(hgamma_eval$fun * fymu_all_eval$d2fun) / mixprob) + w_beta_score^2
  w_gamma2_hess <- ((-(fymu_all_eval$fun - fy) * hgamma_eval$d2fun) / mixprob) + w_gamma_score^2
  w_beta_gamma_hess <- (-(fymu_all_eval$dfun * hgamma_eval$dfun) / mixprob) +
   ((fymu_all_eval$fun - fy) * hgamma_eval$fun * hgamma_eval$dfun * fymu_all_eval$dfun) / (mixprob^2)
  
  Xw4 <- sweep(x, 1, w_beta2_hess, "*")
  Deltaw6 <- sweep(z, 1, w_gamma2_hess, "*")
  Xw5 <- sweep(x, 1, w_beta_gamma_hess, "*")
  
  Hess <- matrix(0, d + ncol(z), d + ncol(z))
  Hess[1:d, 1:d] <- crossprod(x, Xw4)
  Hess[(d+1):ncol(Hess), (d+1):ncol(Hess)] <- crossprod(z, Deltaw6)
  Hess[1:d, (d+1):ncol(Hess)] <- crossprod(Xw5, z)
  Hess[(d+1):ncol(Hess), 1:d] <- t(Hess[1:d, (d+1):ncol(Hess)])
 }
 
 gamma_cur <- as.numeric(-1 * gamma_cur)
 
 cov_1_hat <- tryCatch(solve(Hess, meat), error = function(e) matrix(NA, nrow(Hess), ncol(Hess)))
 covhat <- if (anyNA(cov_1_hat)) matrix(NA, nrow(Hess), ncol(Hess)) else t(solve(Hess, t(cov_1_hat)))
 
 beta_n <- colnames(x)
 if(!is.null(colnames(z))){
  gamma_n <- colnames(z)
 } else{
  gamma_n <- c(1:ncol(z))
 }
 names(gamma_cur) <- gamma_n
 
 if (family_name == "gaussian") {
  if (!anyNA(covhat)) {
   covhat[d+1, ] <- covhat[d+1, ] * (2 * stdcur)
   covhat[, d+1] <- covhat[, d+1] * (2 * stdcur)
  }
  rn <- c(paste("coef", beta_n), "dispersion", paste("m.coef", gamma_n))
 } else {
  rn <- c(paste("coef", beta_n), paste("m.coef", gamma_n))
 }
 rownames(covhat) <- colnames(covhat) <- rn
 
 p_tot <- ncol(covhat)
 p_gamma <- ncol(z)
 idx_other <- 1:(p_tot - p_gamma)
 idx_gamma <- (p_tot - p_gamma + 1):p_tot
 
 if (!anyNA(covhat)) {
  covhat[idx_other, idx_gamma] <- -1 * covhat[idx_other, idx_gamma]
  covhat[idx_gamma, idx_other] <- -1 * covhat[idx_gamma, idx_other]
 }
 
 out <- list(coefficients = beta_cur,
             m.coefficients = gamma_cur,
             residuals = y - mu_cur,
             linear.predictors = eta_cur,
             fitted.values = mu_cur,
             deviance = if (!is.null(w_glm_fit)) w_glm_fit$deviance else NA,
             null.deviance = if (!is.null(w_glm_fit)) w_glm_fit$null.deviance else NA,
             df.residual = n - p,
             df.null = n - 1,
             rank = p,
             family = family,
             converged = converged_flag,
             match.prob = as.numeric(p_cur),
             var = covhat,
             objective = objs[1:iter],
             call = match.call())
 
 if (family_name == "gaussian") out$dispersion <- stdcur^2
 else if (family_name == "gamma") out$dispersion <- 1 / shape_cur
 else out$dispersion <- 1
 
 class(out) <- c("glmMixture")
 return(out)
}

#' Weighted GLM with Mixture-Based Linkage Error Adjustment
#'
#' Fits a generalized linear model (GLM) to linked data with linkage errors,
#' accounting for complex survey designs (weights, stratification, clustering).
#' Implements the pseudo-likelihood approach of Binder (1983) extended to
#' mixture models as described in Slawski et al. (2025).
#'
#' @param formula A symbolic description of the outcome model 
#' (e.g., \code{y ~ x}).
#' @param mismatch.formula A symbolic description of the mismatch model 
#' (e.g., \code{~ z}).
#'   If NULL, an intercept-only model is assumed (constant mismatch rate).
#' @param design A \code{svydesign} object from the \code{survey} package
#'  containing the data, weights, and design information.
#' @param family A description of the error distribution and link function.
#'   Supported: \code{gaussian}, \code{binomial}, \code{poisson}, \code{Gamma}.
#' @param m.rate Optional. A fixed prior for the overall mismatch rate 
#' (proportion 0-1). If provided, it constrains the intercept of the 
#' mismatch model.
#' @param safe.matches Logical vector indicating observations known to be 
#' correct matches.
#' @param control List of control parameters:
#'   \itemize{
#'     \item \code{max.iter}: Max EM iterations (default 1000).
#'     \item \code{tol}: Convergence tolerance (default 1e-6).
#'     \item \code{cmax.iter}: Max iterations for constrained logistic 
#'     (default 1000).
#'   }
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{glmMixtureWeighted} containing coefficients,
#'   robust sandwich standard errors, and design-based variance metrics.
#'
#' @importFrom survey svyvar svymean
#' @importFrom stats glm.fit dnorm dbinom dpois dgamma plogis qlogis density
#'  approxfun
#' @importFrom stats gaussian binomial poisson Gamma weighted.mean
#' @export
glmMixtureWeighted <- function(formula, mismatch.formula = ~ 1, design,
                               family = stats::gaussian, m.rate = NULL,
                               safe.matches = NULL, control = list(), ...) {
 
 con <- list(max.iter = 1000, tol = 1e-6, cmax.iter = 1000)
 con[names(control)] <- control
 
 # Set-up
 data <- design$variables
 # raw weights that sum to N - used for final var Calc (bread & meat)
 raw_weights <- as.vector(stats::weights(design))
 n <- nrow(data)
 
 # optimization weights sum to n: used for glm.fit numerical stability
 weights <- raw_weights * (n / sum(raw_weights))
 
 mf <- stats::model.frame(formula, data = data)
 x <- stats::model.matrix(formula, data = data)
 y <- stats::model.response(mf)
 
 mf_m <- stats::model.frame(mismatch.formula, data = data)
 z <- stats::model.matrix(mismatch.formula, data = data)
 
 if (is.character(family)) {
  if (tolower(family) == "gamma") family <- stats::Gamma(link = "log")
  else family <- get(family, mode = "function", envir = parent.frame())()
 }
 if (is.function(family)) family <- family()
 
 if (is.null(safe.matches)) safe.matches <- rep(FALSE, n)
 
 logitbound <- NULL
 if (!is.null(m.rate)) {
  logitbound <- -log((1 - m.rate) / m.rate)
 }
 
 # Marginal density
 if (family$family == "gaussian") {
  w_mean <- sum(weights * y) / sum(weights)
  w_var  <- sum(weights * (y - w_mean)^2) / (sum(weights) - 1)
  fy <- stats::dnorm(y, mean = w_mean, sd = sqrt(w_var))
 } else if (family$family == "binomial") {
  w_prob <- sum(weights * y) / sum(weights)
  fy <- stats::dbinom(y, size = 1, prob = w_prob)
 } else if (family$family == "poisson") {
  w_lambda <- sum(weights * y) / sum(weights)
  fy <- stats::dpois(y, lambda = w_lambda)
 } else if (family$family == "Gamma") {
  w_mean <- sum(weights * y) / sum(weights)
  w_var <- sum(weights * (y - w_mean)^2) / (sum(weights) - 1)
  shape_marg <- w_mean^2 / w_var
  scale_marg <- w_var / w_mean
  fy <- stats::dgamma(y, shape = shape_marg, scale = scale_marg)
 } else {
  dens <- stats::density(y, weights = weights / sum(weights))
  fy_fun <- stats::approxfun(dens$x, dens$y, rule = 2)
  fy <- fy_fun(y)
 }
 fy[fy < 1e-300] <- 1e-300
 
 # Initialization
 init_fit <- stats::glm.fit(x, y, weights = weights, family = family)
 beta_cur <- stats::coef(init_fit)
 eta <- x %*% beta_cur
 mu <- family$linkinv(eta)
 
 sum_w <- sum(weights)
 p_dim <- ncol(x)
 dispersion_cur <- 1
 shape_cur <- NULL
 
 if (family$family == "gaussian") {
  dispersion_cur <- sum(weights * (y - mu)^2) / (sum_w - p_dim)
 } else if (family$family == "Gamma") {
  dispersion_est <- sum(weights * (y - mu)^2 / mu^2) / (sum_w - p_dim)
  shape_cur <- 1 / dispersion_est
  dispersion_cur <- 1/shape_cur 
 }
 
 gamma_cur <- rep(0, ncol(z))
 if (!is.null(logitbound)) {
  if (ncol(z) == 1 && all(z == 1)) gamma_cur <- rep(-logitbound, ncol(z))
  else gamma_cur <- c(min(-logitbound, 0), rep(0, ncol(z) - 1))
 }
 
 h_func <- function(lp) {
  pi <- stats::plogis(lp)
  list(fun = pi, dfun = pi * (1 - pi), 
       d2fun = pi * (1 - pi) * (1 - 2 * pi))
 }
 
 # EM-Algorithm
 loglik_old <- -Inf
 converged <- FALSE
 history <- numeric(con$max.iter)
 
 for (iter in 1:con$max.iter) {
  
  # E-step
  h_z <- stats::plogis(z %*% gamma_cur)
  
  if (family$family == "gaussian") {
   f_yx <- stats::dnorm(y, mean = mu, sd = sqrt(dispersion_cur))
  } else if (family$family == "binomial") {
   f_yx <- stats::dbinom(y, size = 1, prob = mu)
  } else if (family$family == "poisson") {
   f_yx <- stats::dpois(y, lambda = mu)
  } else if (family$family == "Gamma") {
   f_yx <- stats::dgamma(y, shape = shape_cur, rate = shape_cur / mu)
  }
  
  num <- f_yx * h_z
  denom <- num + fy * (1 - h_z)
  post_prob <- num / denom
  post_prob[safe.matches] <- 1
  post_prob[denom == 0] <- 0
  
  # M-step
  if (!is.null(logitbound)) {
   glm_h <- constrained_weighted_logistic(z, 1 - post_prob, weights, logitbound, con$cmax.iter)
   gamma_cur <- -glm_h$beta
  } else {
   fit_gamma <- stats::glm.fit(z, post_prob, weights = weights,
                               family = stats::quasibinomial(link = "logit"))
   gamma_cur <- stats::coef(fit_gamma)
  }
  
  w_eff_beta <- weights * post_prob
  
  if (family$family == "Gamma") {
   fit_beta <- stats::glm.fit(x, y, weights = w_eff_beta, family = family)
   beta_cur <- stats::coef(fit_beta)
   eta <- x %*% beta_cur
   mu <- family$linkinv(eta)
   
   sum_w_eff <- sum(w_eff_beta)
   term_k <- sum(w_eff_beta * (log(y/mu) - y/mu))
   for (nr_iter in 1:20) {
    sc <- sum_w_eff * (log(shape_cur) + 1 - digamma(shape_cur)) + term_k
    hess <- sum_w_eff * (1/shape_cur - trigamma(shape_cur))
    new_s <- shape_cur - sc/hess
    sf <- 1
    while (new_s <= 1e-4 && sf > 1e-4) { sf <- sf * 0.5; new_s <- shape_cur -
     (sc/hess) * sf }
    if (new_s <= 1e-4) new_s <- 1e-4
    if (abs(new_s - shape_cur) < 1e-6) { shape_cur <- new_s; break }
    shape_cur <- new_s
   }
   dispersion_cur <- 1/shape_cur
  } else {
   fit_beta <- stats::glm.fit(x, y, weights = w_eff_beta, family = family)
   beta_cur <- stats::coef(fit_beta)
   eta <- x %*% beta_cur
   mu <- family$linkinv(eta)
   
   if (family$family == "gaussian") {
    dispersion_cur <- sum(w_eff_beta * (y - mu)^2) / sum(w_eff_beta)
   }
  }
  
  # Check Convergence (Use n-scaled weights)
  h_z_new <- stats::plogis(z %*% gamma_cur)
  if (family$family == "gaussian") {
   f_yx_new <- stats::dnorm(y, mean = mu, sd = sqrt(dispersion_cur))
  } else if (family$family == "binomial") {
   f_yx_new <- stats::dbinom(y, size = 1, prob = mu)
  } else if (family$family == "poisson") {
   f_yx_new <- stats::dpois(y, lambda = mu)
  } else if (family$family == "Gamma") {
   f_yx_new <- stats::dgamma(y, shape = shape_cur, rate = shape_cur / mu)
  }
  
  lik_i <- (h_z_new * f_yx_new) + ((1 - h_z_new) * fy)
  lik_i[safe.matches] <- f_yx_new[safe.matches]
  lik_i[lik_i <= 0] <- 1e-300
  
  loglik_new <- sum(weights * log(lik_i))
  history[iter] <- loglik_new
  
  if (abs(loglik_new - loglik_old) < con$tol * (abs(loglik_old) + 1)) {
   converged <- TRUE
   break
  }
  loglik_old <- loglik_new
 }
 
 # Variance estimation
 h_res <- h_func(z %*% gamma_cur)
 
 get_derivs <- function() {
  mu_prime <- family$mu.eta(eta)
  var_mu <- family$variance(mu)
  
  if (family$family == "gaussian") {
   phi <- stats::dnorm(y, mean = mu, sd = sqrt(dispersion_cur))
   score_eta <- (y - mu) / var_mu * mu_prime 
   dfun_beta <- phi * (y - mu) / dispersion_cur * mu_prime
   
   # Proper second derivative for beta
   d2fun_beta <- phi * ( ( (y - mu)/dispersion_cur *
                            mu_prime )^2 -
                          (mu_prime^2)/dispersion_cur ) 
   
   score_disp <- -1/(2*dispersion_cur) + (y-mu)^2/(2*dispersion_cur^2)
   dfun_disp <- phi * score_disp
   d2fun_disp <- phi * (score_disp^2 + (1/(2*dispersion_cur^2) - (y-mu)^2/(dispersion_cur^3)))
   
   # Compute Gaussian cross-derivative
   term1 <- (-1/(2*dispersion_cur) + (y-mu)^2/(2*dispersion_cur^2)) * ((y-mu)/dispersion_cur)
   term2 <- -(y-mu)/(dispersion_cur^2)
   d2fun_beta_extra <- phi * (term1 + term2) * mu_prime
   
   return(list(fun=phi, dfun_beta=dfun_beta,
               d2fun_beta=d2fun_beta, 
               dfun_extra=dfun_disp, d2fun_extra=d2fun_disp,
               d2fun_beta_extra=d2fun_beta_extra,
               extra_name="dispersion"))
   
  } else if (family$family == "Gamma") {
   phi <- stats::dgamma(y, shape = shape_cur, 
                        rate = shape_cur/mu)
   score_eta <- (shape_cur / mu^2) * (y - mu) * mu_prime
   dfun_beta <- phi * score_eta
   d2fun_beta <- phi * (score_eta^2 - (shape_cur *
                                        mu_prime^2)/mu^2 )
   
   score_shape <- log(shape_cur) + 1 - digamma(shape_cur) +
    log(y/mu) - y/mu
   dfun_shape <- phi * score_shape
   d2fun_shape <- phi * (score_shape^2 + (1/shape_cur - trigamma(shape_cur)))
   
   # Compute Gamma cross-derivative
   d2fun_beta_extra <- phi * score_shape * score_eta + 
    phi * (1/mu^2 * (y - mu) * mu_prime)
   
   return(list(fun=phi, dfun_beta=dfun_beta,
               d2fun_beta=d2fun_beta,
               dfun_extra=dfun_shape,
               d2fun_extra=d2fun_shape,
               d2fun_beta_extra=d2fun_beta_extra,
               extra_name="shape"))
  } else {
   if (family$family == "binomial") phi <- stats::dbinom(y, size=1, prob=mu)
   else phi <- stats::dpois(y, lambda=mu)
   score_eta <- (y - mu) / var_mu * mu_prime
   dfun_beta <- phi * score_eta
   d2fun_beta <- phi * (score_eta^2 - (mu_prime^2)/var_mu)
   
   return(list(fun=phi, dfun_beta=dfun_beta,
               d2fun_beta=d2fun_beta,
               dfun_extra=NULL, 
               d2fun_extra=NULL, d2fun_beta_extra=NULL))
  }
 }
 
 derivs <- get_derivs()
 
 # Meat (J): Calculate unweighted Scores
 denom_i <- h_res$fun * derivs$fun + (1 - h_res$fun) * fy
 denom_i[safe.matches] <- derivs$fun[safe.matches]
 # Safety check
 denom_i[denom_i < 1e-300] <- 1e-300
 
 score_beta_i <- (h_res$fun * derivs$dfun_beta) / denom_i
 U_beta <- sweep(x, 1, as.vector(score_beta_i), "*")
 
 score_gamma_i <- ((derivs$fun - fy) * h_res$dfun) / denom_i
 U_gamma <- sweep(z, 1, as.vector(score_gamma_i), "*")
 
 U_extra <- NULL
 if (!is.null(derivs$dfun_extra)) {
  score_extra_i <- (h_res$fun * derivs$dfun_extra) / denom_i
  U_extra <- matrix(score_extra_i, ncol = 1)
  colnames(U_extra) <- derivs$extra_name
 }
 
 U_total <- cbind(U_beta, U_extra, U_gamma)
 
 # Clean column names for svytotal
 safe_names <- paste0("Sc_", 1:ncol(U_total))
 colnames(U_total) <- safe_names
 
 # Add scores to design to use survey functions
 tmp_des <- design
 tmp_des$variables <- as.data.frame(U_total)
 
 # Calculate J using survey package (handle design weights correctly)
 fmla <- as.formula(paste("~", paste(safe_names, collapse = "+")))
 J <- stats::vcov(survey::svytotal(fmla, tmp_des))
 
 # Bread (H): Calculate using raw weights (population scale)
 # This puts H on the same N-scale as J, ensuring correct cancellation
 t1_bb <- (h_res$fun * derivs$d2fun_beta) / denom_i
 H_scalar_bb <- t1_bb - score_beta_i^2
 
 t1_gg <- ((derivs$fun - fy) * h_res$d2fun) / denom_i
 H_scalar_gg <- t1_gg - score_gamma_i^2
 
 t1_bg <- (derivs$dfun_beta * h_res$dfun) / denom_i
 H_scalar_bg <- t1_bg - (score_gamma_i * score_beta_i)
 
 p <- ncol(x); q <- ncol(z)
 n_extra <- if(is.null(U_extra)) 0 else 1
 dim_h <- p + n_extra + q
 H <- matrix(0, dim_h, dim_h)
 
 idx_b <- 1:p
 idx_e <- if(n_extra) p + 1 else integer(0)
 idx_g <- (p + n_extra + 1):dim_h
 
 W_bb <- as.vector(raw_weights * H_scalar_bb)
 H[idx_b, idx_b] <- crossprod(x, sweep(x, 1, W_bb, "*"))
 
 W_gg <- as.vector(raw_weights * H_scalar_gg)
 H[idx_g, idx_g] <- crossprod(z, sweep(z, 1, W_gg, "*"))
 
 W_bg <- as.vector(raw_weights * H_scalar_bg)
 block_bg <- crossprod(x, sweep(z, 1, W_bg, "*"))
 H[idx_b, idx_g] <- block_bg
 H[idx_g, idx_b] <- t(block_bg)
 
 if (n_extra > 0) {
  t1_ee <- (h_res$fun * derivs$d2fun_extra) / denom_i
  H_scalar_ee <- t1_ee - score_extra_i^2
  H[idx_e, idx_e] <- sum(raw_weights * H_scalar_ee)
  
  t1_be <- (h_res$fun * derivs$d2fun_beta_extra) / denom_i
  H_scalar_be <- t1_be - (score_beta_i * score_extra_i)
  
  block_be <- crossprod(x, as.vector(raw_weights *
                                      H_scalar_be)) 
  H[idx_b, idx_e] <- block_be
  H[idx_e, idx_b] <- t(block_be)
  
  t1_eg <- (derivs$dfun_extra * h_res$dfun) / denom_i
  H_scalar_eg <- t1_eg - (score_extra_i * score_gamma_i)
  block_eg <- crossprod(z, as.vector(raw_weights *
                                      H_scalar_eg)) 
  H[idx_g, idx_e] <- block_eg
  H[idx_e, idx_g] <- t(block_eg)
 }
 
 H_inv <- try(solve(H), silent = TRUE)
 if (inherits(H_inv, "try-error")) H_inv <- MASS::ginv(H)
 
 V_sand <- H_inv %*% J %*% H_inv
 
 all_names <- c(colnames(x), if(n_extra) derivs$extra_name else NULL, 
                colnames(z))
 rownames(V_sand) <- colnames(V_sand) <- all_names
 
 res <- list(coefficients = beta_cur, m.coefficients = gamma_cur,
             dispersion = if(family$family=="gaussian") dispersion_cur else if(family$family=="Gamma") 1/shape_cur else 1,
             vcov = V_sand, se = sqrt(diag(V_sand)),
             naive_vcov = -H_inv * (sum(raw_weights) / n), # naive variance (scale H_inv from N to n)
             naive_se = sqrt(diag(-H_inv * (sum(raw_weights) / n))),
             weights = weights, converged = converged, iterations = iter,
             call = match.call(), family = family)
 class(res) <- "glmMixtureWeighted"
 return(res)
}

constrained_weighted_logistic <- function(X, y, w, bound, cmax.iter) {
 X <- as.matrix(X)
 d <- ncol(X)
 barx <- colSums(w * X) / sum(w)
 nloglik <- function(mu) {
  cutoff <- 1E-10
  mu <- pmax(pmin(mu, 1 - cutoff), cutoff)
  ll <- sum(w * (y * log(mu) + (1 - y) * log(1 - mu)))
  return(-ll)
 }
 grad_Hess <- function(mu) {
  resid <- y - mu
  grad <- -crossprod(X, w * resid)
  h_w <- as.vector(w * mu * (1 - mu))
  Xw <- sweep(X, 1, sqrt(h_w), "*")
  Hess <- crossprod(Xw)
  return(list(grad = grad, Hess = Hess))
 }
 betacur <- c(bound, rep(0, d - 1))
 mu <- stats::plogis(X %*% betacur)
 iter <- 1
 objs <- numeric(cmax.iter)
 objs[1] <- nloglik(mu)
 while (iter < cmax.iter) {
  gH <- grad_Hess(mu)
  rhs <- -gH$grad + (gH$Hess %*% betacur)
  betanew <- try(solve(gH$Hess, rhs), silent=TRUE)
  if(inherits(betanew, "try-error")) betanew <- MASS::ginv(gH$Hess) %*% rhs
  eta_tmp <- X %*% betanew
  mean_eta <- sum(w * eta_tmp) / sum(w)
  if (mean_eta > bound) {
   lhs_aug <- rbind(cbind(gH$Hess, barx), c(barx, 0))
   rhs_aug <- c(rhs, bound)
   sol_aug <- try(solve(lhs_aug, rhs_aug), silent=TRUE)
   if(inherits(sol_aug, "try-error")) sol_aug <- MASS::ginv(lhs_aug) %*% rhs_aug
   betanew <- sol_aug[1:d]
  }
  m <- 0
  upd <- betanew - betacur
  while(m < 20) {
   step_beta <- betacur + upd * (0.5^m)
   mu_step <- stats::plogis(X %*% step_beta)
   obj_step <- nloglik(mu_step)
   if (obj_step < objs[iter] + 1e-8) {
    betacur <- step_beta
    mu <- mu_step
    break
   }
   m <- m + 1
  }
  iter <- iter + 1
  objs[iter] <- nloglik(mu)
  if (abs(objs[iter] - objs[iter-1]) < 1e-8) break
 }
 return(list(beta = betacur))
}

glmMixtureWeighted_Bootstrap <- function(formula, 
                                         mismatch.formula = ~ 1, design,
                                         family = stats::gaussian,
                                         m.rate = NULL,
                                         replicates = 50, ...) {
 
 fit_main <- glmMixtureWeighted(formula, mismatch.formula, design,
                                family, m.rate, ...)
 
 # create replicate design
 # This handles the stratification/clustering automatically
 boot_design <- survey::as.svrepdesign(design, type = "bootstrap", replicates = replicates) 
 
 # function to fit the model on a specific replicate weight set
 run_replicate <- function(w_rep, data_rep) {
  
  # use the original data structure but updated weights
  rep_des <- survey::svydesign(ids = ~1, weights = ~w_rep, 
                               data = data_rep)
  tryCatch({
   fit <- suppressWarnings(glmMixtureWeighted(formula,                                               mismatch.formula, 
                                              rep_des,
                                              family, 
                                              m.rate, ...))
   
   if (any(abs(fit$m.coefficients) > 20, na.rm = TRUE)) return(NULL)
   if (any(abs(fit$coefficients) > 100, na.rm = TRUE)) return(NULL)
   
   return(c(fit$coefficients, fit$dispersion, fit$m.coefficients))
  }, error = function(e) return(NULL))
 }
 
 # Extract weights matrix from replicate design
 rep_weights <- weights(boot_design, type = "analysis")
 
 # Apply function to all replicates
 boot_res_list <- lapply(1:ncol(rep_weights), function(i) run_replicate(rep_weights[,i], boot_design$variables))
 boot_res_list <- boot_res_list[!sapply(boot_res_list, is.null)]
 
 if (length(boot_res_list) < 10) return(fit_main)
 
 # Aggregate results
 boot_res <- do.call(rbind, boot_res_list)
 
 # survey-design bootstrap scaling
 scale_factor <- boot_design$scale
 if (is.null(scale_factor)) scale_factor <- 1
 
 # Centered crossproduct with scale factor attached
 means <- colMeans(boot_res)
 sweep_res <- sweep(boot_res, 2, means, "-")
 cov_boot <- scale_factor * crossprod(sweep_res)
 
 
 fit_main$vcov <- cov_boot
 fit_main$se <- sqrt(diag(cov_boot))
 names(fit_main$se) <- colnames(boot_res)
 return(fit_main)
}

set.seed(2026)

n_trials <- 1000
n_boot <- 50 # Using 50 to balance runtime, can be increased

true_beta <- c(1.0, 2.0, -1.5, 1.0) # Intercept, X1 (Continuous), X2 (Binary), X3 (Interaction)
true_gamma <- c(0.0, 1.0)           # Intercept, Z (Slope of 1.0 mimics match probabilities)

safe_val <- function(x) {
 if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA)
 return(as.numeric(unname(x))[1]) 
}

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id != "") {
 chunk_size <- 20
 end_trial <- as.integer(task_id) * chunk_size
 start_trial <- end_trial - (chunk_size - 1)
 run_trials <- start_trial:end_trial
 results <- NULL
} else {
 run_trials <- 1:n_trials
 results <- NULL
}

for (s in run_trials) {
 print(paste("Running Trial:", s))
 
 set.seed(2026 + s)
 
 # Generate Population: 10 Strata, 100 clusters/stratum, 100 units/cluster (N = 100,000)
 pop_list <- list()
 cluster_id_counter <- 1
 
 for (h in 1:10) {
  for (c in 1:100) {
   n_c <- 100
   
   x1_c <- rnorm(n_c, mean = h * 0.2, sd = 1) 
   x2_c <- rbinom(n_c, 1, 0.5)                
   x3_c <- x1_c * x2_c                        
   
   # Generate covariates for the mismatch model (mimicking match probabilities)
   ps_c <- rbeta(n_c, 4.5, 0.5)
   logit_ps_c <- log(ps_c / (1 - ps_c))
   while(sum(is.infinite(logit_ps_c)) != 0) {
    ps_c <- rbeta(n_c, 4.5, 0.5)
    logit_ps_c <- log(ps_c / (1 - ps_c))
   }
   z_c <- logit_ps_c
   
   u_c <- rnorm(1, 0, 0.25) 
   y_c <- true_beta[1] + true_beta[2]*x1_c + true_beta[3]*x2_c + true_beta[4]*x3_c + u_c + rnorm(n_c, 0, 1)
   
   pop_list[[length(pop_list) + 1]] <- data.frame(
    stratum = h, cluster = cluster_id_counter, unit = 1:n_c,
    y = y_c, x1 = x1_c, x2 = x2_c, x3 = x3_c, z = z_c
   )
   cluster_id_counter <- cluster_id_counter + 1
  }
 }
 pop_df <- do.call(rbind, pop_list)
 
 pop_df$size_meas <- exp(0.2 * pop_df$y + 0.1 * pop_df$x1 + rnorm(nrow(pop_df), 0, 0.1))
 
 samp_list <- list()
 for (h in 1:10) {
  strat_data <- pop_df[pop_df$stratum == h, ]
  
  cluster_sizes <- tapply(strat_data$size_meas, strat_data$cluster, sum)
  pi_c <- 10 * (cluster_sizes / sum(cluster_sizes))
  pi_c <- pmin(pi_c, 1) 
  
  sampled_clusters <- sample(names(cluster_sizes), 10, prob = cluster_sizes)
  
  for (sc in sampled_clusters) {
   cl_data <- strat_data[strat_data$cluster == as.numeric(sc), ]
   pi_u_cond <- 50 * (cl_data$size_meas / sum(cl_data$size_meas))
   pi_u_cond <- pmin(pi_u_cond, 1)
   
   sampled_units <- sample(1:nrow(cl_data), 50, prob = cl_data$size_meas)
   
   samp_data <- cl_data[sampled_units, ]
   samp_data$weights <- 1 / (pi_c[sc] * pi_u_cond[sampled_units])
   samp_list[[length(samp_list) + 1]] <- samp_data
  }
 }
 samp_df <- do.call(rbind, samp_list)
 samp_df$weights_n <- samp_df$weights * (nrow(samp_df) / sum(samp_df$weights))
 
 target_mrates <- c(0.05, 0.15, 0.30)
 
 trial_results <- NULL 
 
 for (rate in target_mrates) {
  
  # find the intercept (gamma_0) for the target mismatch rate
  find_int <- function(int) {
   mean(plogis(int + true_gamma[2] * samp_df$z)) - (1 - rate)
  }
  
  calibrated_intercept <- uniroot(find_int, lower = -10, upper = 10, extendInt = "yes")$root
  current_gamma <- c(calibrated_intercept, true_gamma[2])
  
  # Introduce Mismatch Error using calibrated intercept
  match_eta <- cbind(1, samp_df$z) %*% current_gamma
  match_prob <- plogis(match_eta)
  m <- 1 - rbinom(nrow(samp_df), prob = match_prob, size = 1)
  samp_df$y_linked <- samp_df$y
  if (sum(m) > 1) {
   shuffled_ix <- sample(which(m == 1)) 
   samp_df$y_linked[shuffled_ix] <- samp_df$y_linked[c(shuffled_ix[2:length(shuffled_ix)], 
                                                       shuffled_ix[1])]
  }
  
  # Define Survey Designs & Model Matrices
  des_complex <- svydesign(id = ~cluster, strata = ~stratum, weights = ~weights, data = samp_df)
  des_weights_only <- svydesign(ids = ~1, weights = ~weights, data = samp_df)
  
  X_mat <- model.matrix(~ x1 + x2 + x3, data = samp_df)
  Z_mat <- model.matrix(~ z, data = samp_df) # Z for SNIL
  
  # Fits
  time <- Sys.time()
  fit_oracle <- try(svyglm(y ~ x1 + x2 + x3, design = des_complex), silent=TRUE)
  time_oracle <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_naive_unwt <- try(glm(y_linked ~ x1 + x2 + x3, data = samp_df), silent=TRUE)
  time_naive_unwt <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_naive_wt <- try(glm(y_linked ~ x1 + x2 + x3, data = samp_df, weights = weights_n), silent=TRUE)
  time_naive_wt <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_naive_des <- try(svyglm(y_linked ~ x1 + x2 + x3, design = des_complex), silent=TRUE)
  time_naive_des <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_mix_unwt <- try(glmMixture(x = X_mat, y = samp_df$y_linked, z = Z_mat, family = "gaussian"), silent=TRUE)
  time_mix_unwt <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_mix_wt <- try(glmMixtureWeighted(formula = y_linked ~ x1 + x2 + x3, mismatch.formula = ~ z, design = des_weights_only, family = gaussian), silent=TRUE)
  time_mix_wt <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_mix_des_sand <- try(glmMixtureWeighted(formula = y_linked ~ x1 + x2 + x3, mismatch.formula = ~ z, design = des_complex, family = gaussian), silent=TRUE)
  time_mix_des_sand <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  time <- Sys.time()
  fit_mix_des_boot <- try(glmMixtureWeighted_Bootstrap(formula = y_linked ~ x1 + x2 + x3, mismatch.formula = ~ z, design = des_complex, family = gaussian, replicates = n_boot), silent=TRUE)
  time_mix_des_boot <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  
  # Obtain results 
  extract_stats <- function(fit, method) {
   b_est <- rep(NA_real_, 4); b_se <- rep(NA_real_, 4)
   g_est <- rep(NA_real_, 2); g_se <- rep(NA_real_, 2)
   disp <- NA_real_
   
   p_len <- 4 # Number of beta parameters
   
   if (!inherits(fit, "try-error")) {
    if (method %in% c("oracle", "naive_unwt", "naive_des")) {
     b_est <- coef(fit)[1:4]
     b_se <- sqrt(diag(vcov(fit)))[1:4]
     disp <- summary(fit)$dispersion
    } else if (method == "naive_wt") {
     b_est <- coef(fit)[1:4]
     b_se <- summary(fit)$coefficients[1:4, 2]
     disp <- summary(fit)$dispersion
    } else if (method == "mix_unwt") {
     b_est <- fit$coefficients[1:4]
     disp <- fit$dispersion
     g_est <- fit$m.coefficients[1:2]
     if (!is.null(fit$var)) {
      ses <- sqrt(diag(fit$var))
      b_se <- ses[1:4]
      g_se <- ses[(p_len + 2):(p_len + 3)]
     }
    } else if (method == "mix_wt") {
     b_est <- fit$coefficients[1:4]
     disp <- fit$dispersion
     g_est <- fit$m.coefficients[1:2]
     b_se <- fit$naive_se[1:4]
     g_se <- fit$naive_se[(p_len + 2):(p_len + 3)]
    } else {
     b_est <- fit$coefficients[1:4]
     disp <- fit$dispersion
     g_est <- fit$m.coefficients[1:2]
     b_se <- fit$se[1:4]
     g_se <- fit$se[(p_len + 2):(p_len + 3)]
    }
   }
   
   calc_cov <- function(est, se, true) {
    if(any(is.na(est)) || any(is.na(se))) return(rep(NA_real_, length(est)))
    as.numeric((est - 1.96*se <= true) & (est + 1.96*se >= true))
   }
   
   time_sec <- get(paste0("time_", method))
   
   b_cov <- calc_cov(b_est, b_se, true_beta)
   g_cov <- calc_cov(g_est, g_se, current_gamma)
   
   return(data.frame(
    Target_Rate = rate,
    Method = method,
    B0_Est = safe_val(b_est[1]), B0_Bias = safe_val(b_est[1] - true_beta[1]), B0_SE = safe_val(b_se[1]), B0_Cov = safe_val(b_cov[1]),
    B1_Est = safe_val(b_est[2]), B1_Bias = safe_val(b_est[2] - true_beta[2]), B1_SE = safe_val(b_se[2]), B1_Cov = safe_val(b_cov[2]),
    B2_Est = safe_val(b_est[3]), B2_Bias = safe_val(b_est[3] - true_beta[3]), B2_SE = safe_val(b_se[3]), B2_Cov = safe_val(b_cov[3]),
    B3_Est = safe_val(b_est[4]), B3_Bias = safe_val(b_est[4] - true_beta[4]), B3_SE = safe_val(b_se[4]), B3_Cov = safe_val(b_cov[4]),
    G0_Est = safe_val(g_est[1]), G0_Bias = safe_val(g_est[1] - current_gamma[1]), G0_SE = safe_val(g_se[1]), G0_Cov = safe_val(g_cov[1]),
    G1_Est = safe_val(g_est[2]), G1_Bias = safe_val(g_est[2] - current_gamma[2]), G1_SE = safe_val(g_se[2]), G1_Cov = safe_val(g_cov[2]),
    Dispersion = safe_val(disp),
    Time = time_sec
   ))
  }
  
  res_step <- rbind(extract_stats(fit_oracle, "oracle"),
                    extract_stats(fit_naive_unwt, "naive_unwt"),
                    extract_stats(fit_naive_wt, "naive_wt"),
                    extract_stats(fit_naive_des, "naive_des"),
                    extract_stats(fit_mix_unwt, "mix_unwt"),
                    extract_stats(fit_mix_wt, "mix_wt"),
                    extract_stats(fit_mix_des_sand, "mix_des_sand"),
                    extract_stats(fit_mix_des_boot, "mix_des_boot"))
  
  trial_results <- rbind(trial_results, res_step)
 } 
 
 if (task_id != "") {
  saveRDS(trial_results, file = sprintf("sim_result_task_%04d.rds", s))
 } else {
  results <- rbind(results, trial_results)
 }
}