library(dplyr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(patchwork)
library(kableExtra)

# Load the pre-simulation setup - part1.RData"

# Modify adj_fun_exact to return the Hessian
adj_fun_exact <- function(X, yperm, Delta) {
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
   objective = objs[1:iter],
   hessian = hess
  )
 )
}

set.seed(2026) 
target_rates <- c(0.05, 0.15, 0.30, 0.50)

diagnostic_results <- list()
cor_matrices <- list()
single_trial_results <- list()

param_names <- c("Beta_Int", "Beta_BMI", "Beta_BMI2", "Beta_Age", "Beta_Female", "Beta_BMIAge", 
                 "Sigma", "Gamma_Int", "Gamma_ZFL", "Gamma_Female", "Gamma_Inc")
beta_sig_vars <- param_names[1:7]
gamma_vars <- param_names

for (rate in target_rates) {
 rate_name <- paste0("Rate_", rate)
 smp_rate <- mapped_cohorts[[rate_name]]
 n_sample <- nrow(smp_rate)
 g0_dgp <- gamma_0_calibrated[[rate_name]]
 
 # Generate true mismatch probabilities
 smp_rate <- smp_rate %>%
  mutate(
   logit_pi = g0_dgp + (gamma_Z * Z_FL) + (gamma_female * RIAGENDRFemale) + (gamma_inc * Inc_S),
   pi_i = plogis(logit_pi)
  )
 
 # Introduce simulated mismatches for a single trial
 smp_trial <- smp_rate
 smp_trial$Y_perm <- smp_trial$DXDTOPF
 smp_trial$m_i <- rbinom(n_sample, size = 1, prob = smp_trial$pi_i)
 
 shuffled_ix <- which(smp_trial$m_i == 1)
 if (length(shuffled_ix) > 1) {
  smp_trial$Y_perm[shuffled_ix] <- smp_trial$Y_perm[c(shuffled_ix[-1], shuffled_ix[1])]
 }
 
 rlz_rate <- mean(smp_trial$m_i)
 Y_corr <- smp_trial$Y_perm
 
 # Setup Matrices
 X_mat <- model.matrix(~ BMXBMI + I(BMXBMI^2) + RIDAGEYR + RIAGENDRFemale + BMXBMI:RIDAGEYR, data = smp_trial)
 D_mat <- model.matrix(~ Z_FL + RIAGENDRFemale + Inc_S, data = smp_trial)
 
 # Fit Model
 t0 <- Sys.time()
 mod <- adj_fun_exact(X_mat, Y_corr, D_mat)
 run_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
 
 # Extract and Scale Hessian (googled)
 H_global <- mod$hessian
 H_beta <- mod$hessian[1:ncol(X_mat), 1:ncol(X_mat)]
 
 scale_hessian <- function(H) {
  d <- diag(H)
  d_safe <- pmax(abs(d), 1e-10) 
  D_inv <- diag(1 / sqrt(d_safe))
  return(D_inv %*% H %*% D_inv)
 }
 
 H_global_scaled <- scale_hessian(H_global)
 H_beta_scaled <- scale_hessian(H_beta)
 
 kappa_global <- kappa(H_global_scaled, exact = TRUE)
 kappa_beta <- kappa(H_beta_scaled, exact = TRUE)
 
 # Extract Correlation Matrix
 cov_mat <- solve(mod$hessian)
 cor_mat <- cov2cor(cov_mat)
 rownames(cor_mat) <- param_names
 colnames(cor_mat) <- param_names
 cor_matrices[[rate_name]] <- cor_mat
 
 # Store Diagnostic Results
 diagnostic_results[[length(diagnostic_results) + 1]] <- data.frame(
  Target_Rate = rate,
  Realized_Rate = rlz_rate,
  Kappa_Beta = kappa_beta,
  Kappa_Global = kappa_global,
  Run_Time_Secs = run_time
 )
 
 # Store Parameter Estimates vs Truth
 single_trial_results[[length(single_trial_results) + 1]] <- data.frame(
  Target_Rate = paste0(rate * 100, "%"),
  Parameter = factor(param_names, levels = param_names),
  Estimate = c(mod$betahat, mod$stdhat, mod$gammahat),
  SE = mod$se
 )
}

# Construct the single True Value reference dataframe
true_df <- data.frame(
 Parameter = factor(param_names, levels = param_names),
 Truth = c(
  sprintf("%.3f", as.numeric(coef(oracle_model))), 
  sprintf("%.3f", summary(oracle_model)$sigma),
  "Calibrated*", 
  sprintf("%.3f", c(-gamma_Z, -gamma_female, -gamma_inc))
 )
)

diag_df <- bind_rows(diagnostic_results)

diag_df %>%
 mutate(
  Target_Rate = paste0(Target_Rate * 100, "%"),
  Realized_Rate = sprintf("%.3f", Realized_Rate),
  Kappa_Beta = sprintf("%.2f", Kappa_Beta),
  Kappa_Global = sprintf("%.2f", Kappa_Global)
 ) %>%
 knitr::kable(
  digits = 3,
  col.names = c("Target Rate", "Realized Rate", "Kappa (Beta)", 
                "Kappa (Global)", "Time (s)"),
  caption = "Condition Numbers of the Scaled Observed Information Matrix",
  align = "ccccc"
 )