.libPaths("~/R/hpc_libs")

args <- commandArgs(trailingOnly = TRUE)
setting_name <- args[1]
array_id <- as.integer(args[2])

library(pldamixture) 

# SETTING CONFIGURATIONS
Truebeta <- c(1, 2, -1.5, 1) 
Truesigma <- 0.25 
Truegamma <- c(1, -2.5, 1) 
corr <- TRUE
target_corr <- -0.8 
K_bins <- 1
sigma_e <- 0
link_func <- "logit"
func_form <- "linear"

if (setting_name == "Step1_Null") {
 Truegamma <- c(1, 0, 1); target_corr <- 0
} else if (setting_name == "Step2_PosCorr") {
 target_corr <- 0.8
} else if (setting_name == "Step3_K5") {
 K_bins <- 5
} else if (setting_name == "Step3_K10") {
 K_bins <- 10
} else if (setting_name == "Step3_K20") {
 K_bins <- 20
} else if (setting_name == "Step4_E0.1") {
 sigma_e <- 0.1
} else if (setting_name == "Step4_E0.5") {
 sigma_e <- 0.5
} else if (setting_name == "Step4_E1.0") {
 sigma_e <- 1.0
} else if (setting_name == "Step5A_Probit") {
 link_func <- "probit"
} else if (setting_name == "Step5B_Cloglog") {
 link_func <- "cloglog"
} else if (setting_name == "Step5C_Interact") {
 func_form <- "interact"
}

# EXACT ADJUSTMENT FUNCTIONS 
adj_fun <- function(X, yperm, Delta) {
 n <- nrow(X)
 
 sv_fun <- function(betahat, stdhat, gammahat) {
  muhat <- as.vector(X %*% betahat)
  hs <- as.vector(plogis(Delta %*% gammahat))
  one_minus_hs <- 1 - hs
  sumhs <- sum(one_minus_hs)
  sumhs_safe <- pmax(sumhs, 1E-300)
  w <- tcrossprod(one_minus_hs) / sumhs_safe
  diag(w) <- hs + (one_minus_hs ^ 2) / sumhs_safe
  phi <- dnorm(matrix(yperm, n, n) - matrix(muhat, n, n, byrow = TRUE), 
               sd = stdhat)
  w_phi <- w * phi
  sumphi <- rowSums(w_phi)
  return(
   list(
    muhat = muhat,
    hs = hs,
    sumhs = sumhs,
    stdhat = stdhat,
    w = w,
    phi = phi,
    w_phi = w_phi,
    sumphi = sumphi
   )
  )
 }
 
 obj_fun <- function(thetahat) {
  sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)], thetahat[(ncol(X) + 2):(length(thetahat))])
  return(sum(-log(pmax(sv$sumphi, 1E-300))))
 }
 
 grad_fun <- function(thetahat) {
  betahat <- thetahat[1:ncol(X)]
  stdhat <- thetahat[(ncol(X) + 1)]
  gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
  sv <- sv_fun(betahat, stdhat, gammahat)
  dbeta <- t(-1 / sv$sumphi) %*%
   ((sv$w * sv$phi * (
    matrix(yperm, n, n) -
     matrix(sv$muhat, n, n, byrow = TRUE)
   ) / stdhat ^ 2) %*% X)
  dsigma <- t(-1 / sv$sumphi) %*%
   rowSums(sv$w * sv$phi * ((
    matrix(yperm, n, n) -
     matrix(sv$muhat, n, n, byrow = TRUE)
   ) ^ 2 / stdhat ^ 3 -
    1 / stdhat))
  dw <- numeric(ncol(Delta))
  for (i in 1:ncol(Delta)) {
   dpi <- c(Delta[, i] * exp(Delta %*% gammahat) / 
             (1 + exp(Delta %*% gammahat)) ^2)
   sumdpi <- sum(-dpi)
   dpij <- outer(sv$hs, dpi, "*") + outer(dpi, sv$hs, "*")
   dwij <- (sv$sumhs * (-(outer(dpi, dpi, "+")) + dpij) -
             sumdpi * (1 - (outer(sv$hs, sv$hs, "+")) +
                        outer(sv$hs, sv$hs, "*"))) / sv$sumhs ^ 2
   diag(dwij) <- diag(dwij) + dpi
   dw[i] <- t(-1 / sv$sumphi) %*% rowSums(sv$phi * dwij)
  }
  return(c(dbeta, dsigma, dw))
 }
 
 maxiter <- 1000
 tol <- 1E-4
 iter <- 1
 betahat <- solve(crossprod(X), crossprod(X, yperm))
 gammahat <- rep(0, ncol(Delta))
 stdhat <- sqrt(sum((yperm - X %*% betahat) ^ 2) / (n - ncol(X)))
 
 sv <- sv_fun(betahat, stdhat, gammahat)
 objs <- numeric(maxiter)
 objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))
 
 while (iter < maxiter) {
  lhat <- sv$w_phi / pmax(sv$sumphi, 1E-300)
  
  diag_lhat <- diag(lhat)
  rowSums_lhat <- rowSums(lhat)
  colSums_lhat <- colSums(lhat)
  y_weighted <- as.vector(crossprod(lhat, yperm))
  
  betahat <- as.vector(solve(crossprod(X, colSums_lhat * X), 
                             crossprod(X, y_weighted)))
  muhat <- as.vector(X %*% betahat)
  var_term <- sum(rowSums_lhat * yperm ^ 2) - 2 * sum(muhat * y_weighted) +
   sum(colSums_lhat * muhat ^ 2)
  stdhat <- max(sqrt(var_term / sum(rowSums_lhat)), 1E-4)
  
  fn <- function(gammaval) {
   hs_tmp <- as.vector(plogis(Delta %*% gammaval))
   sumhs_tmp <- sum(1 - hs_tmp)
   w_ii <- hs_tmp + ((1 - hs_tmp) ^ 2) / pmax(sumhs_tmp, 1E-300)
   log_w_ii <- log(pmax(w_ii, 1E-300))
   log_1_minus_hs <- log(pmax(1 - hs_tmp, 1E-300))
   log_sumhs <- log(pmax(sumhs_tmp, 1E-300))
   term1 <- sum(diag_lhat * log_w_ii)
   term2 <- sum((rowSums_lhat - diag_lhat) * log_1_minus_hs)
   term3 <- sum((colSums_lhat - diag_lhat) * log_1_minus_hs)
   term4 <- sum(rowSums_lhat - diag_lhat) * log_sumhs
   return(-(term1 + term2 + term3 - term4))
  }
  gammahat <- optim(gammahat, fn, method = "BFGS")$par
  
  sv <- sv_fun(betahat, stdhat, gammahat)
  iter <- iter + 1
  objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))
  if (objs[iter] + tol > objs[iter - 1])
   break
 }
 
 hess <- optimHess(par = c(betahat, stdhat, gammahat),
                   fn = obj_fun,
                   gr = grad_fun)
 return(
  list(
   betahat = betahat,
   gammahat = gammahat,
   stdhat = stdhat,
   se = sqrt(diag(solve(hess))),
   objective = objs[1:iter]
  )
 )
}


# Exact Oracle Adjustment Function 
adjustment_function <- function(X, yperm, Delta, beta, gamma, sigma,
                                initbgtype, method, louismethod, lmse_version, est_sigma){
 if(missing(louismethod)){
  louismethod <- 0
  lmse_version <- 1
  est_sigma <- 1
 }
 
 y <- yperm # change to y to align w/ existing adjustment approach code
 n <- nrow(X)
 is_flagged <- rep(FALSE,n)
 p <- ncol(X)
 
 maxiter <- 1000
 tol <- 1E-4
 maxiterouter <- 20 
 tolouter <- 1E-4
 
 # INITIALIZE EM ALGORITHM
 
 ## initbeta and initgamma 
 if (initbgtype == "true"){
  betacur <- beta
  gammacur <- gamma
 } else{
  betacur <- solve(crossprod(X), crossprod(X, y))
  gammacur <- rep(0, ncol(Delta))
 }
 
 ## initial fy value
 if (method == 1 | method == 7){
  kde_obj <- density(y)
  g <- approxfun(x = kde_obj$x, y = kde_obj$y, method = "linear")
  fy <- g(y)
 }
 
 if (method == 2){
  g <- function(yval){dnorm(yval, mean = mean(y), sd = sd(y))} # Gaussian density
  fy <- g(y)
 }
 
 if (method == 3){ 
  weights <- 1 - plogis(Delta %*% gamma)
  weights_g <- weights/sum(weights)
  g <- function(yval) sapply(yval, function(z) sum(weights_g*dnorm(z, mean = X%*%beta, sd = sigma)))
  fy <- g(y)
 }
 
 if (method %in% c(4:6)){
  g <- function(yval, betaval, gammaval){
   weights <- 1 - plogis(Delta %*% gammaval)
   weights_g <- weights/sum(weights)
   fy <- sapply(yval, function(z) sum(weights_g*dnorm(z, mean = X%*%betaval, sd = sigma)))
   
   return(list(weights_g = weights_g, fy = fy))
  }
  
  if (method == 4){
   fy <- g(y, beta, gammacur)$fy
   weights_g <- g(y, beta, gammacur)$weights_g
  } 
  
  if (method == 5){
   fy <- g(y, betacur, gammacur)$fy
   weights_g <- g(y, betacur, gammacur)$weights_g
  }
  
  if (method == 6){
   kde_obj <- density(y)
   g_init <- approxfun(x = kde_obj$x, y = kde_obj$y, method = "linear")
   fy <- g_init(y)
  }
  
 }
 
 ## fymu & nloglik 
 fymu <- function(mu, sigma, sub) dnorm((y - mu)[sub], sd = sigma)
 
 fymu_all_normal <- function(mu, std, sub){
  fun <- dnorm((y - mu)[sub], sd = std)
  d_fun_beta <- fun * (y - mu)[sub]/(std^2)
  d_fun_sigma <- fun*((y - mu)[sub]^2/std^3 - 1/std)
  d2_fun_beta <- d_fun_beta*(y - mu)[sub]/(std^2) - fun/(std^2)
  d2_fun_sigma <- d_fun_sigma*((y - mu)[sub]^2/std^3 - 1/std) +
   fun*(1/std^2 - 3*(y - mu)[sub]^2/std^4)
  d2_fun_beta_sigma <- d_fun_sigma*(y - mu)[sub]/(std^2) -
   2*fun*(y - mu)[sub]/(std^3)
  return(list(fun = fun, dfun_beta = d_fun_beta, dfun_sigma = d_fun_sigma,
              d2fun_beta = d2_fun_beta, d2fun_sigma = d2_fun_sigma,
              d2fun_beta_sigma = d2_fun_beta_sigma))
 }
 
 nloglik <- function(mu, sigma, hs){
  sum(-log(hs[!is_flagged] * fymu(mu, sigma, !is_flagged) +
            (1 - hs[!is_flagged])* fy[!is_flagged])) -
   sum(log(fymu(mu, sigma, is_flagged)))
 }
 
 ## hs
 hgamma <- function(eta){
  pi <- plogis(eta)
  return(list(fun = pi, dfun = pi*(1-pi), d2fun = pi*(1-pi)*(1 - 2*pi)))
 }
 
 # RUN EM-ALGORITHM
 if (method %in% c(1:3)){
  iterouter <- 1
  ## hs
  hs <- hgamma(Delta %*% gammacur)$fun
  
  ## stdcur
  if (est_sigma == 1){
   stdcur <- sqrt(sum((y - X %*% betacur)^2)/(n-p))
  } else{
   stdcur <- sigma
  }
  
  ## mucur
  mucur <- X%*%betacur
  
  ## objs
  iter <- 1
  nloglik_cur <- nloglik(mucur, stdcur, hs)
  objs <- numeric(maxiter)
  objs[iter] <- nloglik_cur
  
  ## pcur
  pcur <- rep(0,n)
  while(iter < maxiter){
   num <-  hs[!is_flagged] * fymu(mucur, stdcur,!is_flagged)
   denom <- num + (1-hs[!is_flagged]) * fy[!is_flagged]
   pcur[!is_flagged] <- num/denom
   pcur[is_flagged] <- 1
   
   if(anyNA(pcur)){
    warning("EM algorithm did not converge. NA observation weight(s) occurred in the E-step.")
    
    return(list(betacur = betacur, gammacur = gammacur, objective = objs, 
                avgmmr = 1-mean(hs), maxgweights = NA, 
                pcur = pcur, emcycleiter = iter, standard.errors = "unknown"))
   }
   
   glm_h <- glm(pcur[!is_flagged] ~ Delta[!is_flagged,] - 1,
                family = quasibinomial)
   gammacur <- coef(glm_h)
   
   hs <- hgamma(Delta %*% gammacur)$fun
   
   wglmfit <- glm(y ~ X - 1, family = gaussian, weights = pcur)
   betacur <- coef(wglmfit)
   mucur <- X %*% betacur
   
   if (est_sigma == 1){
    stdcur <- sqrt(weighted.mean((y - mucur)^2, w = pcur))
   }
   
   iter <- iter + 1
   objs[iter] <- nloglik(mucur, stdcur, hs)
   
   if(is.na(objs[iter])){
    warning("EM algorithm did not converge. NA objective value occurred.")
    
    return(list(betacur = betacur, gammacur = gammacur, objective = objs, 
                avgmmr = 1-mean(hs), maxgweights = NA, 
                pcur = pcur, emcycleiter = iter, standard.errors = "unknown"))
   }
   
   if(objs[iter] + tol > objs[iter-1]){
    break
   }    
  }
  objective <- objs
  weights_maxdiff <- NA
  emcycleiter <- iter
 } else{
  iterouter <- 1
  weights_maxdiff <- numeric(maxiterouter)
  emcycleiter <- numeric(maxiterouter)
  while (iterouter < maxiterouter){
   ## hs
   hs <- hgamma(Delta %*% gammacur)$fun
   
   ## stdcur
   if (est_sigma == 1){
    stdcur <- sqrt(sum((y - X %*% betacur)^2)/(n-p))
   } else{
    stdcur <- sigma
   }
   
   ## mucur
   mucur <- X%*%betacur
   
   ## objs
   iter <- 1
   nloglik_cur <- nloglik(mucur, stdcur, hs)
   objs <- numeric(maxiter)
   objs[iter] <- nloglik_cur
   
   ## pcur
   pcur <- rep(0,n)
   while(iter < maxiter){
    num <-  hs[!is_flagged] * fymu(mucur, stdcur,!is_flagged)
    denom <- num + (1-hs[!is_flagged]) * fy[!is_flagged]
    pcur[!is_flagged] <- num/denom
    pcur[is_flagged] <- 1
    
    if(anyNA(pcur)){
     warning("EM algorithm did not converge. NA observation weight(s) occurred in the E-step.")
     emcycleiter[iterouter] <- iter 
     
     if(iterouter == 1){
      objective <- objs
     } else{
      objective <- cbind(objective, objs)
     }
     return(list(betacur = betacur, gammacur = gammacur, objective = objective, 
                 avgmmr = 1-mean(hs), maxgweights = weights_maxdiff[1:(iterouter)], 
                 pcur = pcur, emcycleiter = emcycleiter[1:(iterouter)],
                 standard.errors = "unknown"))
    }
    
    glm_h <- glm(pcur[!is_flagged] ~ Delta[!is_flagged,] - 1,
                 family = quasibinomial)
    gammacur <- coef(glm_h)
    
    hs <- hgamma(Delta %*% gammacur)$fun
    
    wglmfit <- glm(y ~ X - 1, family = gaussian, weights = pcur)
    betacur <- coef(wglmfit)
    mucur <- X %*% betacur
    
    if (est_sigma == 1){
     stdcur <- sqrt(weighted.mean((y - mucur)^2, w = pcur))          
    }
    
    iter <- iter + 1
    objs[iter] <- nloglik(mucur, stdcur, hs)
    
    if(is.na(objs[iter])){
     warning("EM algorithm did not converge. NA objective value occurred.")
     
     if(iterouter == 1){
      objective <- objs
     } else{
      objective <- cbind(objective, objs)
     }
     return(list(betacur = betacur, gammacur = gammacur, objective = objective, 
                 avgmmr = 1-mean(hs), maxgweights = weights_maxdiff[1:(iterouter)], 
                 pcur = pcur, emcycleiter = emcycleiter[1:iterouter],
                 standard.errors = "unknown"))
    }
    
    if(objs[iter] + tol > objs[iter-1]){
     break
    }
   }
   
   # objective update
   if(iterouter == 1){
    objective <- objs
   } else{
    objective <- cbind(objective, objs)
   }
   emcycleiter[iterouter] <- iter
   
   # estimate fy using new coefficient estimate(s)
   if (method == 4){
    newfy <- g(y, beta, gammacur)  
    weights_maxdiff[iterouter] <- max(abs(diff(weights_g - newfy$weights_g)))
    fy <- newfy$fy
    weights_g <- newfy$weights_g
    
    if(weights_maxdiff[iterouter] < tolouter){
     break
    }
   } 
   
   if (method == 5 | (method == 6 & iterouter > 1)){
    newfy <- g(y, betacur, gammacur)  
    weights_maxdiff[iterouter] <- max(abs(diff(weights_g - newfy$weights_g)))
    fy <- newfy$fy
    weights_g <- newfy$weights_g
    
    if(weights_maxdiff[iterouter] < tolouter){
     break
    }
   }
   
   if (method == 6 & iterouter == 1){
    newfy <- g(y, betacur, gammacur)  
    weights_maxdiff[iterouter] <- NA
    fy <- newfy$fy
    weights_g <- newfy$weights_g
   }    
   
   if (method == 7){
    weights <- (1-pcur)
    weights_g <- weights/sum(weights)
    kde_obj <- density(y, weights = weights_g)
    g <- approxfun(x = kde_obj$x, y = kde_obj$y, method = "linear")
    
    fy <- g(y)
    
    if (iterouter == 1){
     weights_maxdiff[iterouter] <- NA
     oldweights <- weights_g
    } else{
     weights_maxdiff[iterouter] <- max(abs(diff(oldweights - weights_g)))
     oldweights <- weights_g
     
     if(weights_maxdiff[iterouter] < tolouter){
      break
     }
    }
   }
   iterouter <- iterouter + 1
  }
 }
 
 # STANDARD ERRORS
 # -------------------------------------------------------------------------
 if (louismethod == 0){
  fymu_all_eval <- fymu_all_normal(mucur, stdcur, !is_flagged)
  
  hgamma_eval <- hgamma(Delta[!is_flagged,] %*% as.matrix(gammacur))
  
  mixprob <- fy[!is_flagged] * (1 - hgamma_eval$fun) +
   hgamma_eval$fun * fymu_all_eval$fun
  
  w_beta_score_num <- (-1) * fymu_all_eval$dfun_beta*hgamma_eval$fun
  w_beta_score_denom <- mixprob
  w_beta_score <- w_beta_score_num/w_beta_score_denom
  
  w_sigma_score_num <- (-1) * fymu_all_eval$dfun_sigma*hgamma_eval$fun
  w_sigma_score_denom <- mixprob
  w_sigma_score <- w_sigma_score_num/w_sigma_score_denom
  
  w_gamma_score_num <- (-1) * (fymu_all_eval$fun -
                                fy[!is_flagged])*hgamma_eval$dfun
  w_gamma_score_denom <- mixprob
  w_gamma_score <- w_gamma_score_num/w_gamma_score_denom
  
  w1 <- w_beta_score^2
  w2 <- w_sigma_score^2
  w3 <- w_gamma_score^2
  
  Xw1 <- sweep(X[!is_flagged,], MARGIN = 1, STATS = w_beta_score, FUN = "*")
  Xw2 <- sweep(matrix(nrow = sum(!is_flagged), ncol = 1, data = 1),
               MARGIN = 1, STATS = w_sigma_score, FUN = "*")
  Deltaw3 <- sweep(as.matrix(Delta[!is_flagged,]), MARGIN = 1,
                   STATS = w_gamma_score, FUN = "*")
  
  meat <- crossprod(cbind(Xw1, Xw2, Deltaw3))
  
  w_beta2_hess <- (-(hgamma_eval$fun * fymu_all_eval$d2fun_beta)/mixprob) + (w_beta_score)^2
  w_sigma2_hess <- (-(hgamma_eval$fun * fymu_all_eval$d2fun_sigma)/mixprob) + (w_sigma_score)^2
  w_gamma2_hess <- ((-(fymu_all_eval$fun - fy[!is_flagged])*hgamma_eval$d2fun)/mixprob) + (w_gamma_score)^2
  w_beta_gamma_hess <- (-(fymu_all_eval$dfun_beta * hgamma_eval$dfun)/mixprob) + ((fymu_all_eval$fun - fy[!is_flagged])*(hgamma_eval$fun)*hgamma_eval$dfun*fymu_all_eval$dfun_beta)/(mixprob^2)
  w_sigma_gamma_hess <- (-(fymu_all_eval$dfun_sigma * hgamma_eval$dfun)/mixprob) + ((fymu_all_eval$fun - fy[!is_flagged])*(hgamma_eval$fun)*hgamma_eval$dfun*fymu_all_eval$dfun_sigma)/(mixprob^2)
  w_beta_sigma_hess <- (-(fymu_all_eval$d2fun_beta_sigma * hgamma_eval$dfun)/mixprob) + w_beta_score * w_sigma_score
  
  Xw4 <- sweep(X[!is_flagged,], MARGIN = 1, STATS = w_beta2_hess, FUN = "*")
  Deltaw6 <- sweep(as.matrix(Delta[!is_flagged,]), MARGIN = 1, STATS = w_gamma2_hess, FUN = "*")
  Xw5 <-  sweep(X[!is_flagged,], MARGIN = 1, STATS = w_beta_gamma_hess, FUN = "*")
  d <- ncol(X)
  Hess <- matrix(nrow = d + 1 + ncol(Delta), ncol =  d + 1 + ncol(Delta))
  one_vector <- matrix(nrow = sum(!is_flagged), ncol = 1, data = 1)
  Hess[1:d, 1:d] <- crossprod(X[!is_flagged,], Xw4)
  Hess[d+1, d+1] <- crossprod(one_vector, sweep(one_vector, MARGIN = 1, STATS = w_sigma2_hess, FUN = "*"))
  Hess[(d+2):(d+ncol(Delta)+1), (d+2):(d+ncol(Delta)+1)] <- crossprod(Delta[!is_flagged,], Deltaw6)
  Hess[1:(d+1), (d+2):(d+ncol(Delta)+1)] <- rbind(crossprod(Xw5, Delta[!is_flagged,]), crossprod(sweep(one_vector, MARGIN = 1, STATS = w_sigma_gamma_hess, FUN = "*"), Delta[!is_flagged,]))
  Hess[(d+2):(d+ncol(Delta)+1), 1:(d+1)] <- t(Hess[1:(d+1), (d+2):(d+ncol(Delta)+1)])
  Hess[1:d, d+1] <- crossprod(X[!is_flagged,], sweep(one_vector, MARGIN = 1, STATS = w_beta_sigma_hess, FUN = "*"))
  Hess[d+1, 1:d] <- t(Hess[1:d, d+1])
  
  cov_1_hat  <- solve(Hess, meat)
  covhat <- t(solve(Hess, t(cov_1_hat)))
  
  ses <- sqrt(diag(covhat))
  nam = NULL
  for (i in 1:(d+ncol(Delta)+1)){
   if(i <= d){nam[i] = paste0("beta", i-1)}
   if(i == d+1){nam[i] = paste0("sigma")}
   if(i>=d+2){nam[i] = paste0("gamma", i-(d+2))}
  }
  names(ses) <- nam
  names(gammacur) <- colnames(Delta)
  
  ses[d+1] <- ses[d+1]*(2*stdcur) # delta method for SE of sigma^2
 } else{
  E_m <- pcur
  Var_m <- pcur * (1 - pcur)
  
  # expected complete data Hessian
  Hc <- matrix(nrow = 4, ncol = 4, data = 0)
  Hc[1:2, 1:2] <- crossprod(sweep(X, MARGIN = 1, FUN = "*", STATS = sqrt(E_m)))/sigma^2
  partial_gamma <- (1-hs) * hs
  Hc[3:4, 3:4] <- crossprod(sweep(Delta, MARGIN = 1, FUN = "*", STATS = sqrt(partial_gamma)))
  
  # covariance of score
  cov_score <- matrix(nrow = 4, ncol = 4, data = 0)
  cov_score[1:2, 1:2] <- crossprod(sweep(X, MARGIN = 1, FUN = "*", STATS = sqrt((y - mucur)^2 * Var_m)))/sigma^4
  
  if(lmse_version == 1){
   cov_score[1:2, 3:4] <- -crossprod(sweep(X, MARGIN = 1, FUN = "*", STATS = E_m * (y - mucur)), sweep(Delta, MARGIN = 1, FUN = "*", STATS = (E_m - hs)))/sigma^2
  } else{
   cov_score[1:2, 3:4] <- -crossprod(sweep(X, MARGIN = 1, FUN = "*", STATS = sqrt(Var_m) * (y - mucur)), sweep(Delta, MARGIN = 1, FUN = "*", STATS = sqrt(Var_m)))/sigma^2
  }
  cov_score[3:4, 1:2] <- t(cov_score[1:2, 3:4])
  cov_score[3:4, 3:4] <- crossprod(sweep(Delta, MARGIN = 1, FUN = "*", STATS = sqrt(Var_m)))
  
  ses <- sqrt(diag(solve(Hc - cov_score)))
  ses <- append(ses, NA, after = length(beta))
 }
 
 return(list(betacur = betacur, gammacur = gammacur, objective = objective, avgmmr = 1-mean(hs), 
             stdcur = stdcur,
             maxgweights = weights_maxdiff[1:(iterouter)], emcycleiter = emcycleiter[1:iterouter], 
             pcur = pcur, standard.errors = ses))
}

# UPDATED DATA GENERATION & MISMATCH ALLOCATION 
data_fun <- function(true_beta, true_gamma, true_sigma, corr, target_corr, 
                     K_bins, sigma_e, link_func, func_form){
 n <- 1000 
 x <- seq(from = -3, to = 3, length = n) 
 d <- c(rep(0, n/2), rep(1, n/2)) 
 X <- cbind(1, d, x, x*d) 
 
 y <- rnorm(n, X %*% true_beta, true_sigma) 
 
 ps <- rbeta(n, 4.5, 0.5)
 logit_ps <- log(ps / (1 - ps))
 while(sum(logit_ps == Inf) != 0){
  ps <- rbeta(n, 4.5, 0.5)
  logit_ps <- log(ps / (1 - ps))
 }
 
 if(corr){
  noise <- rnorm(n)
  int <- target_corr * scale(x) + sqrt(1 - target_corr^2) * noise
  logit_ps <- sort(logit_ps)[rank(int)]
 }
 
 Delta_true <- cbind(1, d, logit_ps)
 Delta_fit <- cbind(1, d, logit_ps + rnorm(n, 0, sigma_e)) # Step 4: Measurement Error
 
 # Step 5: Misspecification of the true linkage model
 if (func_form == "interact") {
  # 5C: Include non-linear interaction in generation, but withhold from Delta_fit
  eta <- (Delta_true %*% true_gamma) + 1.5 * (d * logit_ps)
 } else {
  eta <- Delta_true %*% true_gamma
 }
 
 if (link_func == "probit") {
  prob <- pnorm(eta)
 } else if (link_func == "cloglog") {
  prob <- 1 - exp(-exp(eta))
 } else {
  prob <- plogis(eta) # Standard Logit
 }
 
 m <- 1 - rbinom(n, prob = prob, size = 1) 
 mmr <- sum(m)/n
 
 # Step 3: Stratified Local Shuffling (A1 Violation)
 yperm <- y
 mismatch_idx <- which(m == 1)
 if (length(mismatch_idx) > 1) {
  if (K_bins > 1) {
   x_mis <- x[mismatch_idx]
   sort_order <- order(x_mis)
   sorted_idx <- mismatch_idx[sort_order]
   bin_ids <- rep(1:K_bins, length.out = length(mismatch_idx))
   
   for (k in 1:K_bins) {
    bin_idx <- sorted_idx[bin_ids == k]
    if(length(bin_idx) > 1) { yperm[bin_idx] <- sample(yperm[bin_idx]) }
   }
  } else {
   # Global Shuffle
   shuffled_ix <- sample(mismatch_idx)
   yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
  }
 }
 
 return(list(X = X, y = y, yperm = yperm, Delta = Delta_fit,
             true_beta = true_beta, true_gamma = true_gamma, 
             true_sigma = true_sigma, mmr = mmr))
}

# EXECUTION LOOP
trials <- 25
ncoef <- length(Truebeta) + length(Truegamma) + 1 

# Pre-allocate storage
naive_coef <- matrix(nr=trials, nc=ncoef); oracle_coef <- matrix(nr=trials, nc=ncoef)
ext1_coef <- matrix(nr=trials, nc=ncoef); ext1z_coef <- matrix(nr=trials, nc=ncoef)
adelta_coef <- matrix(nr=trials, nc=ncoef); aint_coef <- matrix(nr=trials, nc=ncoef)
az_coef <- matrix(nr=trials, nc=ncoef); aoracle_coef <- matrix(nr=trials, nc=ncoef)

naive_se <- matrix(nr=trials, nc=ncoef); oracle_se <- matrix(nr=trials, nc=ncoef)
ext1_se <- matrix(nr=trials, nc=ncoef); ext1z_se <- matrix(nr=trials, nc=ncoef)
adelta_se <- matrix(nr=trials, nc=ncoef); aint_se <- matrix(nr=trials, nc=ncoef)
az_se <- matrix(nr=trials, nc=ncoef); aoracle_se <- matrix(nr=trials, nc=ncoef)

# Time storage
naive_time <- matrix(nr=trials, nc=1); oracle_time <- matrix(nr=trials, nc=1)
ext1_time <- matrix(nr=trials, nc=1); ext1z_time <- matrix(nr=trials, nc=1)
adelta_time <- matrix(nr=trials, nc=1); aint_time <- matrix(nr=trials, nc=1)
az_time <- matrix(nr=trials, nc=1); aoracle_time <- matrix(nr=trials, nc=1)

set.seed(array_id)

for (i in 1:trials){
 # REAL-TIME PROGRESS TRACKER
 cat(sprintf("Array %02d | Setting: %s | Running Trial %d of %d\n", array_id, setting_name, i, trials))
 
 data <- data_fun(Truebeta, Truegamma, Truesigma, corr, target_corr, K_bins, sigma_e, link_func, func_form) 
 Y <- data$y; Ystar <- data$yperm; X <- data$X; Delta <- data$Delta; ZDelta <- data$Delta[,-2]
 
 # Naive 
 start <- Sys.time()
 naive <- lm(Ystar ~ X - 1); naive_sigma <- summary(naive)$sigma
 naive_coef[i,] <- c(naive$coefficients, naive_sigma, rep(NA, length(Truegamma)))
 naive_se[i,] <- c(summary(naive)$coefficients[,2], sqrt(2*naive_sigma^4/1000)/(2*naive_sigma), rep(NA, length(Truegamma)))
 stop <- Sys.time()
 naive_time[i,] <- difftime(stop, start, units = "secs")
 
 # Oracle
 start <- Sys.time()
 oracle <- lm(Y ~ X - 1); oracle_sigma <- summary(oracle)$sigma
 oracle_coef[i,] <- c(oracle$coefficients, oracle_sigma, rep(NA, length(Truegamma)))
 oracle_se[i,] <- c(summary(oracle)$coefficients[,2], sqrt(2*oracle_sigma^4/1000)/(2*oracle_sigma), rep(NA, length(Truegamma)))
 stop <- Sys.time()
 oracle_time[i,] <- difftime(stop, start, units = "secs")
 
 # Plain Adjustment (aint)
 start <- Sys.time()
 adj <- fit_mixture(Ystar ~ X - 1) 
 aint_coef[i,] <- c(adj$coefficients, sqrt(adj$dispersion), adj$m.coefficients, rep(NA,2))
 aint_se[i,] <- c(adj$standard.errors[1:(length(Truebeta)+2)], rep(NA, 2))
 stop <- Sys.time()
 aint_time[i,] <- difftime(stop, start, units = "secs")
 
 # Plain Adjustment (adelta)
 start <- Sys.time()
 adjp2 <- fit_mixture(Ystar ~ X - 1, mformula = ~Delta[,-1])
 adelta_coef[i,] <- c(adjp2$coefficients, sqrt(adjp2$dispersion), adjp2$m.coefficients)
 adelta_se[i,] <- adjp2$standard.errors
 stop <- Sys.time()
 adelta_time[i,] <- difftime(stop, start, units = "secs")
 
 # Plain Adjustment (az)
 start <- Sys.time()
 adjp3 <- fit_mixture(Ystar ~ X - 1, mformula = ~ZDelta[,-1])
 az_coef[i,] <- c(adjp3$coefficients, sqrt(adjp3$dispersion), adjp3$m.coefficients[1], NA, adjp3$m.coefficients[2])
 az_se[i,] <- c(adjp3$standard.errors[1:(length(Truebeta)+2)], NA, adjp3$standard.errors[(ncoef-1)])
 stop <- Sys.time()
 az_time[i,] <- difftime(stop, start, units = "secs")
 
 # Extension 1 (ext1)
 start <- Sys.time()
 adje1 <- adj_fun(X, Ystar, Delta = Delta) 
 ext1_coef[i,] <- c(adje1$betahat, adje1$stdhat, adje1$gammahat)
 ext1_se[i,] <- adje1$se 
 stop <- Sys.time()
 ext1_time[i,] <- difftime(stop, start, units = "secs")
 
 # Extension 1 (ext1z)
 start <- Sys.time()
 adje2 <- adj_fun(X, Ystar, Delta = ZDelta) 
 ext1z_coef[i,] <- c(adje2$betahat, adje2$stdhat, adje2$gammahat[1], NA, adje2$gammahat[2])
 ext1z_se[i,] <- c(adje2$se[1:(length(Truebeta)+2)], NA, adje2$se[(ncoef-1)])
 stop <- Sys.time()
 ext1z_time[i,] <- difftime(stop, start, units = "secs")
 
 # Oracle Adjustment (aoracle)
 start <- Sys.time()
 adjo <- try(adjustment_function(X, Ystar, Delta, Truebeta, Truegamma, Truesigma, "naive", 3, 0, 1, 1), silent = TRUE)
 if(!is.character(adjo)){
  aoracle_coef[i,] <- c(adjo$betacur, adjo$stdcur, adjo$gammacur)
  aoracle_se[i,] <- adjo$standard.errors
 }
 stop <- Sys.time()
 aoracle_time[i,] <- difftime(stop, start, units = "secs")
}

# SAVE RESULTS
out_dir <- paste0("results_", setting_name)
if (!dir.exists(out_dir)) dir.create(out_dir)

saveRDS(naive_coef, sprintf("%s/naive_coef-%02d.rds", out_dir, array_id))
saveRDS(oracle_coef, sprintf("%s/oracle_coef-%02d.rds", out_dir, array_id))
saveRDS(aint_coef, sprintf("%s/aint_coef-%02d.rds", out_dir, array_id))
saveRDS(adelta_coef, sprintf("%s/adelta_coef-%02d.rds", out_dir, array_id))
saveRDS(az_coef, sprintf("%s/az_coef-%02d.rds", out_dir, array_id))
saveRDS(ext1_coef, sprintf("%s/ext1_coef-%02d.rds", out_dir, array_id))
saveRDS(ext1z_coef, sprintf("%s/ext1z_coef-%02d.rds", out_dir, array_id))
saveRDS(aoracle_coef, sprintf("%s/aoracle_coef-%02d.rds", out_dir, array_id))

saveRDS(naive_se, sprintf("%s/naive_se-%02d.rds", out_dir, array_id))
saveRDS(oracle_se, sprintf("%s/oracle_se-%02d.rds", out_dir, array_id))
saveRDS(aint_se, sprintf("%s/aint_se-%02d.rds", out_dir, array_id))
saveRDS(adelta_se, sprintf("%s/adelta_se-%02d.rds", out_dir, array_id))
saveRDS(az_se, sprintf("%s/az_se-%02d.rds", out_dir, array_id))
saveRDS(ext1_se, sprintf("%s/ext1_se-%02d.rds", out_dir, array_id))
saveRDS(ext1z_se, sprintf("%s/ext1z_se-%02d.rds", out_dir, array_id))
saveRDS(aoracle_se, sprintf("%s/aoracle_se-%02d.rds", out_dir, array_id))

saveRDS(naive_time, sprintf("%s/naive_time-%02d.rds", out_dir, array_id))
saveRDS(oracle_time, sprintf("%s/oracle_time-%02d.rds", out_dir, array_id))
saveRDS(aint_time, sprintf("%s/aint_time-%02d.rds", out_dir, array_id))
saveRDS(adelta_time, sprintf("%s/adelta_time-%02d.rds", out_dir, array_id))
saveRDS(az_time, sprintf("%s/az_time-%02d.rds", out_dir, array_id))
saveRDS(ext1_time, sprintf("%s/ext1_time-%02d.rds", out_dir, array_id))
saveRDS(ext1z_time, sprintf("%s/ext1z_time-%02d.rds", out_dir, array_id))
saveRDS(aoracle_time, sprintf("%s/aoracle_time-%02d.rds", out_dir, array_id))