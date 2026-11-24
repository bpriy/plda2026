# simulation_exact.R
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.integer(args[1])
setting <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set parameters based on the scenario (1, 2, or 3)
if (scenario == 1) {
 Truegamma <- c(1, -2.5, 1)
 corr <- TRUE
 target_corr <- -0.8
} else if (scenario == 2) {
 Truegamma <- c(1, -2.5, 1)
 corr <- FALSE
 target_corr <- NULL
} else if (scenario == 3) {
 Truegamma <- c(1, 0, 1)
 corr <- TRUE
 target_corr <- -0.8
} else {
 stop("Invalid scenario. Must be 1, 2, or 3.")
}

trials <- 25
Truebeta <- c(1, 2, -1.5, 1)
Truesigma <- 0.25 
ncoef <- length(Truebeta) + 1 + length(Truegamma)

##########################################################
# EXACT ADJUSTMENT FUNCTION
##########################################################
adj_fun_exact <- function(X, yperm, Delta){
 n <- nrow(X)
 
 sv_fun <- function(betahat, stdhat, gammahat){
  muhat <- as.vector(X %*% betahat) 
  hs <- as.vector(plogis(Delta %*% gammahat))
  one_minus_hs <- 1 - hs
  sumhs <- sum(one_minus_hs) 
  sumhs_safe <- pmax(sumhs, 1E-300)                                           
  w <- tcrossprod(one_minus_hs) / sumhs_safe                                      
  diag(w) <- hs + (one_minus_hs^2)/sumhs_safe                                       
  phi <- dnorm(matrix(yperm, n, n) - matrix(muhat, n, n, byrow=TRUE), sd = stdhat)
  w_phi <- w*phi
  sumphi <- rowSums(w_phi) 
  return(list(muhat = muhat, hs = hs, sumhs = sumhs, 
              stdhat = stdhat, w = w, phi = phi, w_phi = w_phi, 
              sumphi = sumphi))
 }
 
 obj_fun <- function(thetahat){
  sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)],
               thetahat[(ncol(X) + 2):(length(thetahat))])
  return(sum(-log(pmax(sv$sumphi, 1E-300)))) 
 }
 
 grad_fun <- function(thetahat){
  betahat <- thetahat[1:ncol(X)]; stdhat <- thetahat[(ncol(X) + 1)]
  gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
  sv <- sv_fun(betahat, stdhat, gammahat)
  dbeta <- t(-1/sv$sumphi) %*% ((sv$w*sv$phi*(matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))/stdhat^2) %*% X)
  dsigma <- t(-1/sv$sumphi) %*% 
   rowSums(sv$w*sv$phi*((matrix(yperm, n, n) - 
                          matrix(sv$muhat, n, n, byrow=TRUE))^2/stdhat^3 -
                         1/stdhat))
  dw <- numeric(ncol(Delta))
  for (i in 1:ncol(Delta)){
   dpi <- c(Delta[,i] * exp(Delta %*% gammahat)/(1+exp(Delta %*% gammahat))^2)
   sumdpi <- sum(-dpi); dpij <- outer(sv$hs, dpi, "*") + outer(dpi, sv$hs, "*")
   dwij <- (sv$sumhs * (-(outer(dpi, dpi, "+")) + dpij) - 
             sumdpi * (1-(outer(sv$hs, sv$hs, "+")) +  
                        outer(sv$hs, sv$hs, "*")))/sv$sumhs^2
   diag(dwij) <- diag(dwij) + dpi
   dw[i] <- t(-1/sv$sumphi) %*% rowSums(sv$phi * dwij)
  }
  return(c(dbeta, dsigma, dw))
 }
 
 maxiter <- 1000; tol <- 1E-4; iter <- 1 
 betahat <- solve(crossprod(X), crossprod(X, yperm))
 gammahat <- rep(0, ncol(Delta))
 stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
 
 sv <- sv_fun(betahat, stdhat, gammahat)
 objs <- numeric(maxiter)
 objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300))) 
 
 while(iter < maxiter){
  lhat <- sv$w_phi / pmax(sv$sumphi, 1E-300)                            
  
  diag_lhat <- diag(lhat)                                                       
  rowSums_lhat <- rowSums(lhat)                                                 
  colSums_lhat <- colSums(lhat)                                                 
  y_weighted <- as.vector(crossprod(lhat, yperm)) 
  
  betahat <- as.vector(solve(crossprod(X, colSums_lhat * X), 
                             crossprod(X, y_weighted)))
  muhat <- as.vector(X %*% betahat) 
  var_term <- sum(rowSums_lhat * yperm^2) - 2 * sum(muhat * y_weighted) +
   sum(colSums_lhat * muhat^2)
  stdhat <- max(sqrt(var_term / sum(rowSums_lhat)), 1E-4) 
  
  fn <- function(gammaval){                                                    
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
  gammahat <- optim(gammahat, fn, method = "BFGS")$par 
  
  sv <- sv_fun(betahat, stdhat, gammahat)
  iter <- iter + 1
  objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))                            
  if(objs[iter] + tol > objs[iter-1]) break
 }
 
 hess <- optimHess(par = c(betahat, stdhat, gammahat), 
                   fn = obj_fun, gr = grad_fun)
 return(list(betahat = betahat, gammahat = gammahat,
             stdhat = stdhat, se = sqrt(diag(solve(hess))),
             objective = objs[1:iter]))
}
##########################################################
# DATA GENERATION FUNCTIONS
##########################################################
mismatches_fun <- function(gammaval, Delta, y, n, p){
 m <- 1 - rbinom(n, prob = plogis(Delta %*% gammaval), size = 1) 
 mmr <- sum(m)/n
 yperm <- y
 shuffled_ix <- sample(which(m == 1)) 
 yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
 return(list(yperm = yperm, mmr = mmr, m = m))
}

data_fun <- function(true_beta, true_gamma, true_sigma, corr = FALSE, target_corr = NULL){
 n <- 1000 
 p <- length(true_beta)
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
 
 Delta <- cbind(matrix(nrow = n, ncol = 1, data = 1), d, logit_ps)
 mismatches <- mismatches_fun(true_gamma, Delta, y, n, p)   
 
 return(list(X = X, y = y, yperm = mismatches$yperm, Delta = Delta))
}

##########################################################
# SIMULATION LOOP
##########################################################
aexact_coef <- matrix(nr = trials, nc = ncoef, data = NA)
aexact_se   <- matrix(nr = trials, nc = ncoef, data = NA)
aexact_time <- matrix(nr = trials, nc = 1, data = NA)

set.seed(setting)

for (i in 1:trials){
 print(i)
 data <- data_fun(true_beta = Truebeta, true_gamma = Truegamma, 
                  true_sigma = Truesigma, corr = corr, target_corr = target_corr) 
 
 start <- Sys.time()
 res <- try(adj_fun_exact(X = data$X, yperm = data$yperm, Delta = data$Delta))
 stop <- Sys.time()
 
 if(inherits(res, "try-error")) next
 
 aexact_coef[i, ] <- c(res$betahat, res$stdhat, res$gammahat)
 aexact_se[i, ]   <- res$se
 aexact_time[i, ] <- difftime(stop, start, units = "secs")
}

##########################################################
# SAVE OUTPUT
##########################################################
out_dir <- sprintf("results_s411_exact_scenario%d", scenario)
if (!dir.exists(out_dir)) dir.create(out_dir)

saveRDS(aexact_coef, sprintf("%s/aexact_coef-%02d.rds", out_dir, setting))
saveRDS(aexact_se,   sprintf("%s/aexact_se-%02d.rds", out_dir, setting))
saveRDS(aexact_time, sprintf("%s/aexact_time-%02d.rds", out_dir, setting))