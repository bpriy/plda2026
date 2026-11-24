.libPaths("~/R/hpc_libs")

library(Matrix)
library(stats)

# HPC Setup & Parallelization Logic
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) task_id <- 1

batch_id <- task_id
trials_per_batch <- 10
set.seed(2026 + batch_id)

cat(sprintf("Starting Large-N Scaling Task: %d | Batch: %d | Trials: %d\n", 
            task_id, batch_id, trials_per_batch))

# 5. Exact O(n) + KDE
adj_fun_exactkde <- function(X, yperm, Delta){
 n <- nrow(X)
 fast_kde <- function(eval_points, centers, bw, weights, fixed_bounds = NULL) {
  if(all(weights == 0)) return(rep(0, length(eval_points)))
  if(any(weights < 0)) {
   w_pos <- pmax(weights, 0)
   w_neg <- pmax(-weights, 0)
   res_pos <- fast_kde(eval_points, centers, bw, w_pos, fixed_bounds)
   res_neg <- fast_kde(eval_points, centers, bw, w_neg, fixed_bounds)
   return(res_pos - res_neg)
  }
  sum_w <- sum(weights)
  if (sum_w == 0) return(rep(0, length(eval_points)))
  if(is.null(fixed_bounds)) {
   min_val <- min(c(centers, eval_points)) - 4 * bw
   max_val <- max(c(centers, eval_points)) + 4 * bw
  } else {
   min_val <- fixed_bounds[1]
   max_val <- fixed_bounds[2]
  }
  dens <- stats::density(centers, bw = bw, weights = weights / sum_w, 
                         n = 2048, from = min_val, to = max_val, 
                         kernel = "gaussian")
  approx_vals <- stats::splinefun(dens$x, dens$y, method = "natural")(eval_points)
  
  return(approx_vals * sum_w)
 }
 obj_fun <- function(thetahat, fixed_bounds = NULL){
  betahat_tmp <- thetahat[1:ncol(X)]
  stdhat_tmp <- thetahat[(ncol(X) + 1)]
  gammahat_tmp <- thetahat[(ncol(X) + 2):(length(thetahat))]
  
  muhat_tmp <- as.vector(X %*% betahat_tmp) 
  hs_tmp <- as.vector(plogis(Delta %*% gammahat_tmp))
  v_tmp <- 1 - hs_tmp
  S_safe <- pmax(sum(v_tmp), 1E-300)
  kde1_vals <- fast_kde(eval_points = yperm, centers = muhat_tmp, 
                        bw = stdhat_tmp, weights = v_tmp, 
                        fixed_bounds = fixed_bounds)
  diag_phi <- dnorm(yperm - muhat_tmp, sd = stdhat_tmp)
  sumphi_tmp <- (v_tmp / S_safe) * kde1_vals + hs_tmp * diag_phi
  
  return(sum(-log(pmax(sumphi_tmp, 1E-300)))) 
 }
 
 maxiter <- 1000; tol <- 1E-4; iter <- 1 
 betahat <- as.vector(solve(crossprod(X), crossprod(X, yperm)))
 gammahat <- rep(0, ncol(Delta))
 stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
 
 objs <- numeric(maxiter)
 objs[iter] <- obj_fun(c(betahat, stdhat, gammahat))
 
 while(iter < maxiter){
  muhat <- as.vector(X %*% betahat) 
  hs <- as.vector(plogis(Delta %*% gammahat))
  v <- 1 - hs
  S_safe <- pmax(sum(v), 1E-300)
  diag_phi <- dnorm(yperm - muhat, sd = stdhat)
  kde1 <- fast_kde(eval_points = yperm, centers = muhat, 
                   bw = stdhat, weights = v)
  sumphi <- (v / S_safe) * kde1 + hs * diag_phi
  sumphi_safe <- pmax(sumphi, 1E-300)
  w_ii <- hs + (v^2) / S_safe
  diag_lhat <- w_ii * diag_phi / sumphi_safe
  rowSums_lhat <- rep(1, n)
  a_weights <- v / sumphi_safe
  kde2 <- fast_kde(eval_points = muhat, centers = yperm, 
                   bw = stdhat, weights = a_weights)
  colSums_lhat <- (v / S_safe) * kde2 + (hs * diag_phi) / sumphi_safe
  b_weights <- (v * yperm) / sumphi_safe
  kde3 <- fast_kde(eval_points = muhat, centers = yperm, 
                   bw = stdhat, weights = b_weights)
  y_weighted <- (v / S_safe) * kde3 + ((hs * diag_phi) / sumphi_safe) * yperm
  betahat <- as.vector(solve(crossprod(X, colSums_lhat * X), 
                             crossprod(X, y_weighted)))
  muhat_new <- as.vector(X %*% betahat) 
  var_term <- sum(rowSums_lhat * yperm^2) - 2 * sum(muhat_new * y_weighted) + 
   sum(colSums_lhat * muhat_new^2)
  stdhat <- max(sqrt(var_term / sum(rowSums_lhat)), 1E-4) 
  
  fn_gamma <- function(gammaval){                                                    
   hs_tmp <- as.vector(plogis(Delta %*% gammaval))                             
   sumhs_tmp <- sum(1 - hs_tmp)                                                
   w_ii <- hs_tmp + ((1 - hs_tmp)^2) / pmax(sumhs_tmp, 1E-300)                 
   log_w_ii <- log(pmax(w_ii, 1E-300))                                         
   log_1_minus_hs <- log(pmax(1 - hs_tmp, 1E-300))                             
   log_sumhs <- log(pmax(sumhs_tmp, 1E-300))                                   
   term1 <- sum(diag_lhat * log_w_ii)                                          
   term2 <- sum((rowSums_lhat - diag_lhat) * log_1_minus_hs)                   
   term3 <- sum((colSums_lhat - diag_lhat) * log_1_minus_hs)                   
   term4 <- sum(rowSums_lhat - diag_lhat) * log_sumhs                          
   return(-(term1 + term2 + term3 - term4))                                    
  } 
  gammahat <- optim(gammahat, fn_gamma, method = "BFGS")$par 
  iter <- iter + 1
  objs[iter] <- obj_fun(c(betahat, stdhat, gammahat))                            
  if(objs[iter] + tol > objs[iter-1]) break
 }
 final_muhat <- as.vector(X %*% betahat)
 fixed_grid_bounds <- c(min(c(final_muhat, yperm)) - 4 * stdhat, 
                        max(c(final_muhat, yperm)) + 4 * stdhat)
 
 obj_fun_fixed <- function(thetahat) {
  obj_fun(thetahat, fixed_bounds = fixed_grid_bounds)
 }
 hess <- optimHess(par = c(betahat, stdhat, gammahat), fn = obj_fun_fixed)
 
 return(list(betahat = betahat, gammahat = gammahat,
             stdhat = stdhat, se = sqrt(diag(solve(hess))),
             objective = objs[1:iter]))
}

# Helper
check_monotonic <- function(obj) { all(diff(obj) <= 1e-8) }

################################################################################
# Simulation Logic
Truebeta <- c(1, -1)
Truegamma <- c(2.5, 4.5) 
Truesigma <- 0.25

N_seq <- c(1000, 10000, 100000, 1000000)
results_scaling <- list()

for (r in 1:trials_per_batch) {
 
 n_list <- list()
 
 for (n in N_seq) {
  cat(sprintf("[%s] Running Trial %d of %d | N = %d\n", 
              Sys.time(), r, trials_per_batch, n))
  flush.console()
  
  # 1. Data Generation
  X1 <- rnorm(n)
  X <- model.matrix(~ X1)
  Delta <- X
  L <- c(X %*% Truebeta)
  Y <- rnorm(n, L, Truesigma)
  
  # 2. Block Logic
  n_blocks <- 10
  blocks_vector <- rep(1:n_blocks, each = n / n_blocks)
  m <- 1 - rbinom(n, prob = plogis(Delta %*% Truegamma), size = 1)
  Ystar <- Y
  
  # 3. Shuffle exactly within blocks
  for(b in 1:n_blocks){
   shuffled_ix <- sample(x = which(m == 1 & blocks_vector == b), 
                         size = sum(m == 1 & blocks_vector == b))
   if(length(shuffled_ix) >= 2) Ystar[shuffled_ix] <- Ystar[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
  }
  
  # 4. Estimation
  t1 <- Sys.time()
  res <- adj_fun_exactkde(X, Ystar, Delta)
  rt <- as.numeric(difftime(Sys.time(), t1, units="secs"))
  
  n_list[[as.character(n)]] <- data.frame(
   Batch = batch_id, Trial = r, N = n, Method = "Exact + KDE",
   Beta_0 = res$betahat[1], Beta_1 = res$betahat[2],
   SE_Beta_0 = res$se[1], SE_Beta_1 = res$se[2],
   Gamma_0 = res$gammahat[1], Gamma_1 = res$gammahat[2],
   SE_Gamma_0 = res$se[4], SE_Gamma_1 = res$se[5],
   Sigma = res$stdhat, SE_Sigma = res$se[3],
   Runtime_Secs = rt, Iters = length(res$objective) - 1,
   Monotonic = check_monotonic(res$objective)
  )
 }
 
 results_scaling[[r]] <- do.call(rbind, n_list)
}

final_results <- do.call(rbind, results_scaling)
out_dir <- "results_scalability/N_scaling"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_path <- sprintf("%s/batch_results_%03d.rds", out_dir, batch_id)
saveRDS(final_results, save_path)
cat(sprintf("Saved Large-N Scaling batch results to %s\n", save_path))