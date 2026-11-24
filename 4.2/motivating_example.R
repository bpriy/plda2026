set.seed(42)

# 1. Define Population Parameters
N <- 100000
s <- 0.30        # Swap rate
alpha <- 0.20    # Mismatch rate

# Generate Confounder Z (50/50 split)
Z <- sample(1:2, N, replace = TRUE)
U <- runif(N)
X <- rep(0, N)
Y <- rep(0, N)

# Generate X and Y (Correlated, but not perfectly, dependent on Z)
# Z = 1: High prevalence. P(1,1)=0.6, P(1,0)=0.1, P(0,1)=0.1, P(0,0)=0.2
idx_z1 <- which(Z == 1)
X[idx_z1] <- ifelse(U[idx_z1] >= 0.3, 1, 0)
Y[idx_z1] <- ifelse(U[idx_z1] >= 0.4 | (U[idx_z1] >= 0.2 & U[idx_z1] < 0.3), 1, 0)

# Z = 2: Low prevalence. P(1,1)=0.2, P(1,0)=0.1, P(0,1)=0.1, P(0,0)=0.6
idx_z2 <- which(Z == 2)
X[idx_z2] <- ifelse(U[idx_z2] >= 0.7, 1, 0)
Y[idx_z2] <- ifelse(U[idx_z2] >= 0.8 | (U[idx_z2] >= 0.6 & U[idx_z2] < 0.7), 1, 0)

# 2. Oracle Estimator (P_true)
P_oracle <- prop.table(table(X, Y))
print(round(P_oracle, 3))

# 3. Stage 1: SDC (Data Swapping within Z)
Y_tilde <- Y
for (z in 1:2) {
 idx <- which(Z == z)
 n_z <- length(idx)
 swap_idx <- sample(idx, size = round(s * n_z))
 # Swap with a random draw from the exact empirical margin of Y within Z
 Y_tilde[swap_idx] <- sample(Y[idx], size = length(swap_idx), replace = TRUE)
}

# 4. Stage 2: Record Linkage (Mismatch Error)
Y_star <- Y_tilde
mismatch_idx <- sample(1:N, size = round(alpha * N))
# Mismatches drawn from global empirical margin of Y_tilde
Y_star[mismatch_idx] <- sample(Y_tilde, size = length(mismatch_idx), replace = TRUE)

# 5. Naive Estimator (P_obs)
C <- table(X, Y_star) # Observed counts
P_naive <- prop.table(C)
print(round(P_naive, 3))

# 6. Pre-compute Fixed Noise Matrices (P_ind and P_CI)
p_x <- prop.table(table(X))
p_y_star <- prop.table(table(Y_star))
P_ind <- outer(p_x, p_y_star)

P_CI <- matrix(0, nrow = 2, ncol = 2)
for (z in 1:2) {
 p_z <- sum(Z == z) / N
 p_x_z <- prop.table(table(X[Z == z]))
 p_y_star_z <- prop.table(table(Y_star[Z == z]))
 P_CI <- P_CI + outer(p_x_z, p_y_star_z) * p_z
}

#############################################################
# Estimator 1: Plain Adjustment EM (Matching Slawski et al.)
#############################################################
p_sat <- P_naive
# Initialization: Independence component based on start point margins
p_ind_plain <- outer(rowSums(p_sat), colSums(p_sat))
eps <- .Machine$double.eps

calc_nloglik <- function(p_s, p_i, rate) {
 -sum(C * log(p_s * (1 - rate) + p_i * rate + eps))
}

objs_plain <- numeric(1000)
objs_plain[1] <- calc_nloglik(p_sat, p_ind_plain, alpha)

for(iter in 2:1000) {
 # E-Step
 numer <- p_sat * (1 - alpha)
 denom <- numer + (p_ind_plain * alpha)
 weights <- numer / (denom + eps)

 # M-Step: P_sat
 counts_sat <- C * weights
 sum_sat <- sum(counts_sat)
 if (sum_sat > 0) p_sat <- counts_sat / sum_sat

 # M-Step: P_ind
 counts_ind <- C * (1 - weights)
 sum_ind <- sum(counts_ind)
 if (sum_ind > 1e-12) {
  p_ind_plain <- outer(rowSums(counts_ind) / sum_ind, colSums(counts_ind) / sum_ind)
 }

 # Objective Tracking & Convergence
 objs_plain[iter] <- calc_nloglik(p_sat, p_ind_plain, alpha)
 if (abs(objs_plain[iter] - objs_plain[iter - 1]) < 1e-6) break
}

theta_plain <- p_sat
print(paste("Plain Adjustment EM converged in", iter, "iterations ---"))
print(round(theta_plain, 3))
all(diff(objs_plain[1:iter]) <= 0) # Verification of monotonic decrease

#############################################################
# Estimator 2: Plain Adjustment with Gamma
#############################################################
gamma <- alpha + s - (alpha * s)
p_sat_gamma <- P_naive
p_ind_gamma <- outer(rowSums(p_sat_gamma), colSums(p_sat_gamma))

objs_gamma <- numeric(1000)
objs_gamma[1] <- calc_nloglik(p_sat_gamma, p_ind_gamma, gamma)

for (iter in 2:1000) {
 # E-Step
 numer <- p_sat_gamma * (1 - gamma)
 denom <- numer + (p_ind_gamma * gamma)
 weights <- numer / (denom + eps)

 # M-Step
 counts_sat <- C * weights
 sum_sat <- sum(counts_sat)
 if (sum_sat > 0) p_sat_gamma <- counts_sat / sum_sat

 counts_ind <- C * (1 - weights)
 sum_ind <- sum(counts_ind)
 if (sum_ind > 1e-12) {
  p_ind_gamma <- outer(rowSums(counts_ind) / sum_ind, colSums(counts_ind) / sum_ind)
 }

 # Objective Tracking & Convergence
 objs_gamma[iter] <- calc_nloglik(p_sat_gamma, p_ind_gamma, gamma)
 if (abs(objs_gamma[iter] - objs_gamma[iter - 1]) < 1e-6) break
}

print(paste("Plain Adjustment (Gamma) EM converged in", iter, "iterations ---"))
print(round(p_sat_gamma, 3))
all(diff(objs_gamma[1:iter]) <= 0) # Verification of monotonic decrease

#############################################################
# Estimator 3: Proposed Three-Component EM
#############################################################
# Fixed weights defined by the swapping and mismatch mechanisms
w1 <- (1 - alpha) * (1 - s)
w2 <- (1 - alpha) * s
w3 <- alpha

theta_prop <- P_naive

# NLL specific to the 3-component model
calc_nloglik_prop <- function(th) {
 -sum(C * log(w1 * th + w2 * P_CI + w3 * P_ind + eps))
}

objs_prop <- numeric(1000)
objs_prop[1] <- calc_nloglik_prop(theta_prop)

for (iter in 2:1000) {
 # E-Step
 denominator <- (w1 * theta_prop) + (w2 * P_CI) + (w3 * P_ind) + eps
 tau_prop <- (w1 * theta_prop) / denominator

 # M-Step
 # Note: P_CI and P_ind are fixed, so no update required here
 C_tilde_prop <- C * tau_prop
 theta_prop <- C_tilde_prop / sum(C_tilde_prop)

 # Objective Tracking & Convergence
 objs_prop[iter] <- calc_nloglik_prop(theta_prop)
 if (abs(objs_prop[iter] - objs_prop[iter - 1]) < 1e-6) break
}

print(paste("Proposed Mixture EM converged in", iter, "iterations ---"))
print(round(theta_prop, 3))
all(diff(objs_prop[1:iter]) <= 0) # Verification of monotonic decrease
