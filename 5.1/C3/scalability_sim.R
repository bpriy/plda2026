.libPaths("~/R/hpc_libs")

library(Matrix)
library(stats)

# HPC Setup & Parallelization Logic
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) task_id <- 1

# Split the 300 array tasks into 3 distinct parallel jobs
if (task_id <= 100) {
  job_type <- "methods"
  batch_id <- task_id
} else if (task_id <= 200) {
  job_type <- "correct"
  batch_id <- task_id - 100
} else {
  job_type <- "misspec"
  batch_id <- task_id - 200
}

trials_per_batch <- 10
set.seed(2026 + batch_id)

cat(sprintf("Starting Task: %d | Job Type: %s | Batch: %d | Trials: %d\n", 
            task_id, job_type, batch_id, trials_per_batch))

# 6 Method Implementations 
# 1. Original Dense O(n^2)
adj_fun_original <- function(X, yperm, Delta){
  n <- nrow(X)
  sv_fun <- function(betahat, stdhat, gammahat){
    muhat <- as.vector(X %*% betahat) 
    hs <- as.vector(plogis(Delta %*% gammahat))
    one_minus_hs <- 1 - hs
    sumhs <- sum(one_minus_hs) 
    sumhs_safe <- pmax(sumhs, 1E-300)
    w <- tcrossprod(one_minus_hs) / sumhs_safe
    diag(w) <- hs + (one_minus_hs^2)/sumhs_safe
    diff_mat <- matrix(yperm, n, n) - matrix(muhat, n, n, byrow=TRUE)
    phi <- dnorm(diff_mat, sd = stdhat)
    w_phi <- w*phi
    sumphi <- rowSums(w_phi) 
    return(list(muhat = muhat, hs = hs, sumhs = sumhs, stdhat = stdhat, w = w, phi = phi, w_phi = w_phi, sumphi = sumphi))
  }
  obj_fun <- function(thetahat){
    sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)], thetahat[(ncol(X) + 2):(length(thetahat))])
    return(sum(-log(pmax(sv$sumphi, 1E-300))))
  }
  grad_fun <- function(thetahat){
    betahat <- thetahat[1:ncol(X)]; stdhat <- thetahat[(ncol(X) + 1)]; gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
    sv <- sv_fun(betahat, stdhat, gammahat)
    dbeta <- t(-1/sv$sumphi) %*% ((sv$w*sv$phi*(matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))/stdhat^2) %*% X)
    dsigma <- t(-1/sv$sumphi) %*% rowSums(sv$w*sv$phi*((matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))^2/stdhat^3 - 1/stdhat))
    dw <- numeric(ncol(Delta))
    for (i in 1:ncol(Delta)){
      dpi <- c(Delta[,i] * exp(Delta %*% gammahat)/(1+exp(Delta %*% gammahat))^2)
      sumdpi <- sum(-dpi); dpij <- outer(sv$hs, dpi, "*") + outer(dpi, sv$hs, "*")
      dwij <- (sv$sumhs * (-(outer(dpi, dpi, "+")) + dpij) - sumdpi * (1-(outer(sv$hs, sv$hs, "+")) + outer(sv$hs, sv$hs, "*")))/sv$sumhs^2
      diag(dwij) <- diag(dwij) + dpi
      dw[i] <- t(-1/sv$sumphi) %*% rowSums(sv$phi * dwij)
    }
    return(c(dbeta, dsigma, dw))
  }
  maxiter <- 1000; tol <- 1E-4; iter <- 1 
  betahat <- solve(crossprod(X), crossprod(X, yperm))
  gammahat <- rep(0, ncol(Delta)); stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
  sv <- sv_fun(betahat, stdhat, gammahat)
  objs <- numeric(maxiter); objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300))) 
  while(iter < maxiter){
    lhat <- sv$w_phi / pmax(sv$sumphi, 1E-300)
    R_i <- rowSums(lhat); C_j <- colSums(lhat); y_weighted <- as.vector(crossprod(lhat, yperm)) 
    betahat <- as.vector(solve(crossprod(X, C_j * X), crossprod(X, y_weighted)))
    muhat <- as.vector(X %*% betahat) 
    var_term <- sum(R_i * yperm^2) - 2 * sum(muhat * y_weighted) + sum(C_j * muhat^2)
    stdhat <- max(sqrt(var_term / sum(R_i)), 1E-4) 
    fn <- function(gammaval){
      hs_tmp <- as.vector(plogis(Delta %*% gammaval))
      one_minus_hs <- 1 - hs_tmp
      sumhs_safe <- pmax(sum(one_minus_hs), 1E-300)
      w_tmp <- tcrossprod(one_minus_hs) / sumhs_safe
      diag(w_tmp) <- hs_tmp + (one_minus_hs^2)/sumhs_safe
      return(-sum(lhat * log(pmax(w_tmp, 1E-300))))
    } 
    gammahat <- optim(gammahat, fn, method = "BFGS")$par 
    sv <- sv_fun(betahat, stdhat, gammahat); iter <- iter + 1
    objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))
    if(objs[iter] + tol > objs[iter-1]) break
  }
  hess <- optimHess(par = c(betahat, stdhat, gammahat), fn = obj_fun, gr = grad_fun)
  return(list(betahat = betahat, gammahat = gammahat, stdhat = stdhat, se = sqrt(diag(solve(hess))), objective = objs[1:iter]))
}

# 2. Sparse (Thresholding)
adj_fun_sparse <- function(X, yperm, Delta){
  n <- nrow(X)
  sv_fun <- function(betahat, stdhat, gammahat){
    muhat <- as.vector(X %*% betahat); hs <- as.vector(plogis(Delta %*% gammahat))
    one_minus_hs <- 1 - hs; sumhs <- sum(one_minus_hs); sumhs_safe <- pmax(sumhs, 1E-300)
    w <- tcrossprod(one_minus_hs) / sumhs_safe; diag(w) <- hs + (one_minus_hs^2)/sumhs_safe
    diff_mat <- matrix(yperm, n, n) - matrix(muhat, n, n, byrow=TRUE); phi <- dnorm(diff_mat, sd = stdhat)
    w_phi <- w*phi; sumphi <- rowSums(w_phi) 
    return(list(muhat=muhat, hs=hs, sumhs=sumhs, stdhat=stdhat, w=w, phi=phi, w_phi=w_phi, sumphi=sumphi))
  }
  obj_fun <- function(thetahat){
    sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)], thetahat[(ncol(X) + 2):(length(thetahat))])
    return(sum(-log(pmax(sv$sumphi, 1E-300))))
  }
  grad_fun <- function(thetahat){
    betahat <- thetahat[1:ncol(X)]; stdhat <- thetahat[(ncol(X) + 1)]; gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
    sv <- sv_fun(betahat, stdhat, gammahat)
    dbeta <- t(-1/sv$sumphi) %*% ((sv$w*sv$phi*(matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))/stdhat^2) %*% X)
    dsigma <- t(-1/sv$sumphi) %*% rowSums(sv$w*sv$phi*((matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))^2/stdhat^3 - 1/stdhat))
    dw <- numeric(ncol(Delta))
    for (i in 1:ncol(Delta)){
      dpi <- c(Delta[,i] * exp(Delta %*% gammahat)/(1+exp(Delta %*% gammahat))^2)
      sumdpi <- sum(-dpi); dpij <- outer(sv$hs, dpi, "*") + outer(dpi, sv$hs, "*")
      dwij <- (sv$sumhs * (-(outer(dpi, dpi, "+")) + dpij) - sumdpi * (1-(outer(sv$hs, sv$hs, "+")) + outer(sv$hs, sv$hs, "*")))/sv$sumhs^2
      diag(dwij) <- diag(dwij) + dpi
      dw[i] <- t(-1/sv$sumphi) %*% rowSums(sv$phi * dwij)
    }
    return(c(dbeta, dsigma, dw))
  }
  maxiter <- 1000; tol <- 1E-4; iter <- 1 
  betahat <- solve(crossprod(X), crossprod(X, yperm))
  gammahat <- rep(0, ncol(Delta)); stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
  sv <- sv_fun(betahat, stdhat, gammahat)
  objs <- numeric(maxiter); objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300))) 
  while(iter < maxiter){
    lhat <- sv$w_phi / pmax(sv$sumphi, 1E-300)
    lhat[lhat <= 1E-5] <- 0 
    R_i <- rowSums(lhat); C_j <- colSums(lhat); y_weighted <- as.vector(crossprod(lhat, yperm)) 
    betahat <- as.vector(solve(crossprod(X, C_j * X), crossprod(X, y_weighted)))
    muhat <- as.vector(X %*% betahat) 
    var_term <- sum(R_i * yperm^2) - 2 * sum(muhat * y_weighted) + sum(C_j * muhat^2)
    stdhat <- max(sqrt(var_term / sum(R_i)), 1E-4) 
    sel_matrix <- lhat > 0
    fn <- function(gammaval){
      hs_tmp <- as.vector(plogis(Delta %*% gammaval))
      one_minus_hs <- 1 - hs_tmp; sumhs_safe <- pmax(sum(one_minus_hs), 1E-300)
      w_tmp <- tcrossprod(one_minus_hs) / sumhs_safe
      diag(w_tmp) <- hs_tmp + (one_minus_hs^2)/sumhs_safe
      return(-sum(lhat[sel_matrix] * log(pmax(w_tmp[sel_matrix], 1E-300)))) 
    } 
    gammahat <- optim(gammahat, fn, method = "BFGS")$par 
    sv <- sv_fun(betahat, stdhat, gammahat); iter <- iter + 1
    objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))
    if(objs[iter] + tol > objs[iter-1]) break
  }
  hess <- optimHess(par = c(betahat, stdhat, gammahat), fn = obj_fun, gr = grad_fun)
  return(list(betahat = betahat, gammahat = gammahat, stdhat = stdhat, se = sqrt(diag(solve(hess))), objective = objs[1:iter]))
}

# 3. Block-Diagonal
adj_fun_block <- function(X, yperm, Delta, blocks){ 
 n <- nrow(X)
 
 blocks_factor <- as.factor(blocks)
 blocks_int <- as.integer(blocks_factor)
 A <- Matrix::fac2sparse(blocks_factor) 
 
 BlockMask <- Matrix::crossprod(A)
 BlockMask <- as(as(BlockMask, "generalMatrix"), "TsparseMatrix") 
 idx_i <- BlockMask@i + 1L
 idx_j <- BlockMask@j + 1L
 diag_mask <- (idx_i == idx_j)
 
 sv_fun <- function(betahat, stdhat, gammahat){
  muhat <- as.numeric(X %*% betahat)
  hs <- as.numeric(plogis(Delta %*% gammahat))
  v <- 1 - hs
  
  S_vec <- as.numeric(Matrix::crossprod(A, A %*% v))
  S_vec <- pmax(S_vec, 1E-300)
  
  diff_val <- yperm[idx_i] - muhat[idx_j]
  phi_val <- dnorm(diff_val, sd = stdhat) 
  
  w_val <- v[idx_i] * v[idx_j] / S_vec[idx_i]
  w_val[diag_mask] <- hs[idx_i[diag_mask]] + 
   (v[idx_i[diag_mask]]^2) / S_vec[idx_i[diag_mask]]
  
  w_phi_val <- w_val * phi_val
  
  L_mat <- Matrix::sparseMatrix(i = idx_i, j = idx_j, 
                                x = w_phi_val, dims = c(n, n))
  sumphi <- as.numeric(Matrix::rowSums(L_mat))
  
  list(sumphi = sumphi, w_phi_val = w_phi_val, phi_val = phi_val, 
       hs = hs, v = v, S_vec = S_vec)
 }
 
 obj_fun <- function(thetahat){ 
  sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)], 
               thetahat[(ncol(X) + 2):(length(thetahat))])
  sum(-log(pmax(sv$sumphi, 1E-300))) 
 }
 
 grad_fun <- function(thetahat){
  betahat <- thetahat[1:ncol(X)]; stdhat <- thetahat[(ncol(X) + 1)]
  gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
  sv <- sv_fun(betahat, stdhat, gammahat)
  
  inv_sumphi <- -1 / pmax(sv$sumphi, 1E-300)
  inv_sumphi_i <- inv_sumphi[idx_i]
  
  muhat <- as.numeric(X %*% betahat)
  diff_val <- yperm[idx_i] - muhat[idx_j]
  
  inner_val_beta <- sv$w_phi_val * diff_val / stdhat^2
  M_beta <- Matrix::sparseMatrix(i = idx_i, j = idx_j, 
                                 x = inner_val_beta, dims = c(n, n))
  dbeta <- as.numeric(Matrix::crossprod(inv_sumphi, M_beta %*% X))
  
  dsigma <- sum(inv_sumphi_i * sv$w_phi_val * 
                 (diff_val^2 / stdhat^3 - 1 / stdhat))
  
  dw <- numeric(ncol(Delta))
  for (k in seq_len(ncol(Delta))){
   dpi <- Delta[,k] * sv$hs * sv$v
   sumdpi_vec <- as.numeric(Matrix::crossprod(A, A %*% (-dpi)))
   
   dpi_i <- dpi[idx_i]; dpi_j <- dpi[idx_j]
   hs_i <- sv$hs[idx_i]; hs_j <- sv$hs[idx_j]
   
   dpij <- hs_i * dpi_j + dpi_i * hs_j
   dwij_val <- ( sv$S_vec[idx_i] * (-(dpi_i + dpi_j) + dpij) - 
                  sumdpi_vec[idx_i] * (1 - (hs_i + hs_j) + hs_i * hs_j) ) / 
    (sv$S_vec[idx_i]^2)
   dwij_val[diag_mask] <- dwij_val[diag_mask] + dpi_i[diag_mask]
   
   dw[k] <- sum(inv_sumphi_i * sv$phi_val * dwij_val)
  }
  c(dbeta, dsigma, dw)
 }
 
 betahat <- as.numeric(solve(crossprod(X), crossprod(X, yperm)))
 gammahat <- rep(0, ncol(Delta))
 stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
 sv <- sv_fun(betahat, stdhat, gammahat)
 iter <- 1; maxiter <- 1000; tol <- 1E-4 
 
 objs <- numeric(maxiter)
 objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))
 
 while(iter < maxiter){
  lhat_val <- sv$w_phi_val / pmax(sv$sumphi[idx_i], 1E-300)
  L_mat <- Matrix::sparseMatrix(i = idx_i, j = idx_j, x = lhat_val, dims = c(n, n))
  
  C_j <- as.numeric(Matrix::colSums(L_mat))
  R_i <- as.numeric(Matrix::rowSums(L_mat))
  y_weighted <- as.numeric(Matrix::crossprod(L_mat, yperm))
  
  diag_lhat_full <- numeric(n)
  diag_lhat_full[idx_i[diag_mask]] <- lhat_val[diag_mask]
  
  betahat <- as.numeric(solve(crossprod(X, C_j * X), crossprod(X, y_weighted)))
  muhat_new <- as.numeric(X %*% betahat) 
  stdhat <- max(sqrt((sum(R_i * yperm^2) - 2 * sum(muhat_new * y_weighted) + 
                       sum(C_j * muhat_new^2)) / sum(R_i)), 1E-4) 
  
  term4_agg <- as.numeric(A %*% (R_i - diag_lhat_full))
  
  fn <- function(gammaval){
   hs_tmp <- as.numeric(plogis(Delta %*% gammaval))
   sumhs_tmp <- as.numeric(A %*% (1 - hs_tmp))
   S_b_mapped <- pmax(sumhs_tmp[blocks_int], 1E-300)
   w_ii <- hs_tmp + ((1 - hs_tmp)^2) / S_b_mapped
   
   term1 <- sum(diag_lhat_full * log(pmax(w_ii, 1E-300)))
   term23 <- sum((R_i + C_j - 2 * diag_lhat_full) * log(pmax(1 - hs_tmp, 1E-300)))
   term4 <- sum(term4_agg * log(pmax(sumhs_tmp, 1E-300)))
   return(-(term1 + term23 - term4))                                   
  }
  gammahat <- optim(gammahat, fn, method = "BFGS")$par 
  
  sv_new <- sv_fun(betahat, stdhat, gammahat)
  iter <- iter + 1; objs[iter] <- sum(-log(pmax(sv_new$sumphi, 1E-300)))
  if(objs[iter] + tol > objs[iter-1]) break
  sv <- sv_new
 }
 
 hess <- optimHess(par = c(betahat, stdhat, gammahat), fn = obj_fun, gr = grad_fun)
 covhat <- solve(hess) 
 
 list(betahat = betahat, gammahat = gammahat, stdhat = stdhat, se = sqrt(diag(covhat)), 
      objective = objs[1:iter])
}

# 4. Exact O(n) Analytical
adj_fun_exact <- function(X, yperm, Delta){
  n <- nrow(X)
  sv_fun <- function(betahat, stdhat, gammahat){
    muhat <- as.vector(X %*% betahat); hs <- as.vector(plogis(Delta %*% gammahat))
    one_minus_hs <- 1 - hs; sumhs <- sum(one_minus_hs); sumhs_safe <- pmax(sumhs, 1E-300)                                           
    w <- tcrossprod(one_minus_hs) / sumhs_safe; diag(w) <- hs + (one_minus_hs^2)/sumhs_safe                                       
    phi <- dnorm(matrix(yperm, n, n) - matrix(muhat, n, n, byrow=TRUE), sd = stdhat)
    w_phi <- w*phi; sumphi <- rowSums(w_phi) 
    return(list(muhat=muhat, hs=hs, sumhs=sumhs, stdhat=stdhat, w=w, phi=phi, w_phi=w_phi, sumphi=sumphi))
  }
  obj_fun <- function(thetahat){
    sv <- sv_fun(thetahat[1:ncol(X)], thetahat[(ncol(X) + 1)], thetahat[(ncol(X) + 2):(length(thetahat))])
    return(sum(-log(pmax(sv$sumphi, 1E-300)))) 
  }
  grad_fun <- function(thetahat){
    betahat <- thetahat[1:ncol(X)]; stdhat <- thetahat[(ncol(X) + 1)]; gammahat <- thetahat[(ncol(X) + 2):(length(thetahat))]
    sv <- sv_fun(betahat, stdhat, gammahat)
    dbeta <- t(-1/sv$sumphi) %*% ((sv$w*sv$phi*(matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))/stdhat^2) %*% X)
    dsigma <- t(-1/sv$sumphi) %*% rowSums(sv$w*sv$phi*((matrix(yperm, n, n) - matrix(sv$muhat, n, n, byrow=TRUE))^2/stdhat^3 - 1/stdhat))
    dw <- numeric(ncol(Delta))
    for (i in 1:ncol(Delta)){
      dpi <- c(Delta[,i] * exp(Delta %*% gammahat)/(1+exp(Delta %*% gammahat))^2)
      sumdpi <- sum(-dpi); dpij <- outer(sv$hs, dpi, "*") + outer(dpi, sv$hs, "*")
      dwij <- (sv$sumhs * (-(outer(dpi, dpi, "+")) + dpij) - sumdpi * (1-(outer(sv$hs, sv$hs, "+")) + outer(sv$hs, sv$hs, "*")))/sv$sumhs^2
      diag(dwij) <- diag(dwij) + dpi
      dw[i] <- t(-1/sv$sumphi) %*% rowSums(sv$phi * dwij)
    }
    return(c(dbeta, dsigma, dw))
  }
  
  maxiter <- 1000; tol <- 1E-4; iter <- 1 
  betahat <- solve(crossprod(X), crossprod(X, yperm))
  gammahat <- rep(0, ncol(Delta)); stdhat <- sqrt(sum((yperm - X %*% betahat)^2)/(n-ncol(X)))
  sv <- sv_fun(betahat, stdhat, gammahat)
  objs <- numeric(maxiter); objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300))) 
  
  while(iter < maxiter){
    lhat <- sv$w_phi / pmax(sv$sumphi, 1E-300)                            
    diag_lhat <- diag(lhat); rowSums_lhat <- rowSums(lhat); colSums_lhat <- colSums(lhat)                                                 
    y_weighted <- as.vector(crossprod(lhat, yperm)) 
    
    betahat <- as.vector(solve(crossprod(X, colSums_lhat * X), crossprod(X, y_weighted)))
    muhat <- as.vector(X %*% betahat) 
    var_term <- sum(rowSums_lhat * yperm^2) - 2 * sum(muhat * y_weighted) + sum(colSums_lhat * muhat^2)
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
    sv <- sv_fun(betahat, stdhat, gammahat); iter <- iter + 1
    objs[iter] <- sum(-log(pmax(sv$sumphi, 1E-300)))                            
    if(objs[iter] + tol > objs[iter-1]) break
  }
  hess <- optimHess(par = c(betahat, stdhat, gammahat), fn = obj_fun, gr = grad_fun)
  return(list(betahat=betahat, gammahat=gammahat, stdhat=stdhat, se=sqrt(diag(solve(hess))), objective=objs[1:iter]))
}

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

# Helper to check monotonicity
check_monotonic <- function(obj) { all(diff(obj) <= 1e-8) }

################################################################################
# Simulation Logic

n <- 1000 
Truebeta <- c(1, -1)
Truegamma <- c(2.5, 4.5) 
Truesigma <- 0.25

results_methods <- list()
results_correct <- list()
results_misspec <- list()

# PART 1: Compare 6 Implementations (10 blocks)
if (job_type == "methods") {
  
  for (r in 1:trials_per_batch) {
  cat(sprintf("[%s] Running Trial %d of %d...\n", Sys.time(), r, trials_per_batch))
    flush.console()
    X1 <- rnorm(n); X <- model.matrix(~ X1); Delta <- X
    L <- c(X %*% Truebeta); Y <- rnorm(n, L, Truesigma)
    
    n_blocks <- 10
    blocks_vector <- rep(1:n_blocks, each = n / n_blocks)
    m <- 1 - rbinom(n, prob = plogis(Delta %*% Truegamma), size = 1)
    Ystar <- Y
    for(b in 1:n_blocks){
      shuffled_ix <- sample(x = which(m == 1 & blocks_vector == b), size = sum(m == 1 & blocks_vector == b))
      if(length(shuffled_ix) >= 2) Ystar[shuffled_ix] <- Ystar[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
    }
    
    run_method <- function(name, func, ...) {
      t1 <- Sys.time(); res <- func(...)
      rt <- as.numeric(difftime(Sys.time(), t1, units="secs"))
      data.frame(
        Batch=batch_id, Trial=r, Method=name,
        Beta_0 = res$betahat[1], Beta_1 = res$betahat[2],
        SE_Beta_0 = res$se[1], SE_Beta_1 = res$se[2],
        Gamma_0 = res$gammahat[1], Gamma_1 = res$gammahat[2],
        SE_Gamma_0 = res$se[4], SE_Gamma_1 = res$se[5],
        Sigma = res$stdhat, SE_Sigma = res$se[3],
        Runtime = rt, Iters = length(res$objective) - 1,
        Monotonic = check_monotonic(res$objective)
      )
    }
    
    m1 <- run_method("Original Dense O(n^2)", adj_fun_original, X, Ystar, Delta)
    m2 <- run_method("Sparse (Thresholding)", adj_fun_sparse, X, Ystar, Delta)
    m3 <- run_method("Block-Diagonal", adj_fun_block, X, Ystar, Delta, blocks=blocks_vector)
    m4 <- run_method("Exact O(n) Analytical", adj_fun_exact, X, Ystar, Delta)
    m5 <- run_method("Exact + KDE", adj_fun_exactkde, X, Ystar, Delta)
    
    results_methods[[r]] <- do.call(rbind, list(m1, m2, m3, m4, m5))
  }
  
  final_results <- do.call(rbind, results_methods)
  out_dir <- "results_scalability/methods"
  
# PART 2: Correctly Specified Blocks 
} else if (job_type == "correct") {

  blocks_seq <- c(1, 10, 20, 50, 100)
  for (r in 1:trials_per_batch) {
  cat(sprintf("[%s] Running Trial %d of %d...\n", Sys.time(), r, trials_per_batch))
    flush.console()
    cr_list <- list()
    for (B in blocks_seq) {
      X1 <- rnorm(n); X <- model.matrix(~ X1); Delta <- X
      L <- c(X %*% Truebeta); Y <- rnorm(n, L, Truesigma)
      blocks_vector <- rep(1:B, each = n / B)
      m <- 1 - rbinom(n, prob = plogis(Delta %*% Truegamma), size = 1)
      Ystar <- Y
      for(b in 1:B){
        shuffled_ix <- sample(x = which(m == 1 & blocks_vector == b), size = sum(m == 1 & blocks_vector == b))
        if(length(shuffled_ix) >= 2) Ystar[shuffled_ix] <- Ystar[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
      }
      
      t_start <- Sys.time()
      res <- adj_fun_block(X, Ystar, Delta, blocks = blocks_vector)
      rt <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
      
      cr_list[[as.character(B)]] <- data.frame(
        Batch = batch_id, Trial = r, Blocks = B, Pairs_Evaluated = (n/B)^2 * B,
        Beta_0 = res$betahat[1], Beta_1 = res$betahat[2],
        SE_Beta_0 = res$se[1], SE_Beta_1 = res$se[2],
        Gamma_0 = res$gammahat[1], Gamma_1 = res$gammahat[2],
        SE_Gamma_0 = res$se[4], SE_Gamma_1 = res$se[5],
        Sigma = res$stdhat, SE_Sigma = res$se[3],
        Runtime_Secs = rt, Iterations = length(res$objective) - 1
      )
    }
    results_correct[[r]] <- do.call(rbind, cr_list)
  }
  
  final_results <- do.call(rbind, results_correct)
  out_dir <- "results_scalability/correct"

# PART 3: Misspecified Blocks
} else if (job_type == "misspec") {

  blocks_fit_seq <- c(1, 10, 20, 50, 100)
  for (r in 1:trials_per_batch) {
  cat(sprintf("[%s] Running Trial %d of %d...\n", Sys.time(), r, trials_per_batch))
    flush.console()
    X1 <- rnorm(n); X <- model.matrix(~ X1); Delta <- X
    L <- c(X %*% Truebeta); Y <- rnorm(n, L, Truesigma)
    true_blocks <- 10
    blocks_true <- rep(1:true_blocks, each = n / true_blocks)
    m <- 1 - rbinom(n, prob = plogis(Delta %*% Truegamma), size = 1)
    Ystar <- Y
    for(b in 1:true_blocks){
      shuffled_ix <- sample(x = which(m == 1 & blocks_true == b), size = sum(m == 1 & blocks_true == b))
      if(length(shuffled_ix) >= 2) Ystar[shuffled_ix] <- Ystar[c(shuffled_ix[2:length(shuffled_ix)], shuffled_ix[1])]
    }
    
    ms_list <- list()
    for (B_fit in blocks_fit_seq) {
      blocks_artificial <- rep(1:B_fit, each = n / B_fit)
      t_start <- Sys.time()
      res <- adj_fun_block(X, Ystar, Delta, blocks = blocks_artificial)
      rt <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
      
      ms_list[[as.character(B_fit)]] <- data.frame(
        Batch = batch_id, Trial = r, Assumed_Blocks = B_fit, Pairs_Evaluated = (n/B_fit)^2 * B_fit,
        Beta_0 = res$betahat[1], Beta_1 = res$betahat[2],
        SE_Beta_0 = res$se[1], SE_Beta_1 = res$se[2],
        Gamma_0 = res$gammahat[1], Gamma_1 = res$gammahat[2],
        SE_Gamma_0 = res$se[4], SE_Gamma_1 = res$se[5],
        Sigma = res$stdhat, SE_Sigma = res$se[3],
        Time_Secs = rt, Iter = length(res$objective) - 1
      )
    }
    results_misspec[[r]] <- do.call(rbind, ms_list)
  }
  
  final_results <- do.call(rbind, results_misspec)
  out_dir <- "results_scalability/misspec"
}

# Save the batch
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_path <- sprintf("%s/batch_results_%03d.rds", out_dir, batch_id)
saveRDS(final_results, save_path)
cat(sprintf("Saved %s batch results to %s\n", job_type, save_path))