# EXACT DATA GENERATION FUNCTION 
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
 Delta_fit <- cbind(1, d, logit_ps + rnorm(n, 0, sigma_e)) 
 
 if (func_form == "interact") {
  eta <- (Delta_true %*% true_gamma) + 1.5 * (d * logit_ps)
 } else {
  eta <- Delta_true %*% true_gamma
 }
 
 if (link_func == "probit") {
  prob <- pnorm(eta)
 } else if (link_func == "cloglog") {
  prob <- 1 - exp(-exp(eta))
 } else {
  prob <- plogis(eta) 
 }
 
 m <- 1 - rbinom(n, prob = prob, size = 1) 
 mmr <- sum(m)/n
 
 return(mmr)
}

settings_to_run <- c("Step1_Null", "Step2_PosCorr")

# Storage for results
mmr_results <- data.frame(Setting = character(), 
                          Average_MMR = numeric(), 
                          stringsAsFactors = FALSE)

for (setting_name in settings_to_run) {
 
 # Base Configurations
 Truebeta <- c(1, 2, -1.5, 1) 
 Truesigma <- 0.25 
 Truegamma <- c(1, -2.5, 1) 
 corr <- TRUE
 target_corr <- -0.8 
 K_bins <- 1
 sigma_e <- 0
 link_func <- "logit"
 func_form <- "linear"
 
 # Apply Setting-Specific Overrides
 if (setting_name == "Step1_Null") {
  Truegamma <- c(1, 0, 1)
  target_corr <- 0
 } else if (setting_name == "Step2_PosCorr") {
  target_corr <- 0.8
 }
 
 # Mimick HPC Array Structure
 num_arrays <- 40
 trials_per_array <- 25
 total_trials <- num_arrays * trials_per_array
 
 mmr_vector <- numeric(total_trials)
 counter <- 1
 
 for (array_id in 1:num_arrays) {
  # Set seed identically to the HPC array environment
  set.seed(array_id)
  
  for (i in 1:trials_per_array) {
   mmr_vector[counter] <- data_fun(Truebeta, Truegamma, Truesigma, corr, 
                                   target_corr, K_bins, sigma_e, 
                                   link_func, func_form)
   counter <- counter + 1
  }
 }
 
 # Store the average mismatch rate
 avg_mmr <- mean(mmr_vector)
 mmr_results <- rbind(mmr_results, data.frame(Setting = setting_name, 
                                              Average_MMR = avg_mmr))
}

print(mmr_results)