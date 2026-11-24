library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)

set.seed(2026)
################################################################################
# Data Generation & Perturbation (Section 4.1)
################################################################################
# Generates the latent true population based on a log-linear model
generate_population <- function(N = 5000) {
 # Confounder Z (e.g., Geographic stratum)
 P_Z <- c(0.3, 0.4, 0.3)
 Z <- sample(1:3, size = N, replace = TRUE, prob = P_Z)

 X <- numeric(N)
 Y <- numeric(N)

 # Log-linear structure ensuring confounding by Z and positive correlation
 for (z in 1:3) {
  idx <- which(Z == z)
  n_z <- length(idx)

  # Base probabilities varying by stratum to induce confounding
  if (z == 1) { p_x <- 0.7; p_y <- 0.7; cor_str <- 0.8 }
  if (z == 2) { p_x <- 0.5; p_y <- 0.5; cor_str <- 0.9 }
  if (z == 3) { p_x <- 0.3; p_y <- 0.3; cor_str <- 0.7 }

  U1 <- runif(n_z)
  U2 <- runif(n_z)

  X_z <- ifelse(U1 < p_x, 1, 0)
  Y_z <- ifelse(runif(n_z) < cor_str, X_z, ifelse(U2 < p_y, 1, 0))

  X[idx] <- X_z
  Y[idx] <- Y_z
 }
 return(data.frame(Z = Z, X = X, Y = Y))
}

# Applies SDC (Stratified Data Swapping) and Linkage Mismatch
perturb_data <- function(data, s, alpha) {
 N <- nrow(data)
 data_obs <- data

 # Stage 1: Disclosure Limitation (M_swap) within Z
 for (z in unique(data$Z)) {
  idx <- which(data$Z == z)
  n_z <- length(idx)
  swap_size <- round(s * n_z)

  if (swap_size > 0) {
   swap_idx <- sample(idx, size = swap_size, replace = FALSE)
   # Draw from exact conditional empirical margin
   data_obs$Y[swap_idx] <- sample(data$Y[idx], size = swap_size, replace = TRUE)
  }
 }

 # Stage 2: Linkage Error (M_link)
 mismatch_size <- round(alpha * N)
 if (mismatch_size > 0) {
  mismatch_idx <- sample(1:N, size = mismatch_size, replace = FALSE)
  # Mismatches drawn from global empirical margin
  data_obs$Y[mismatch_idx] <- sample(data_obs$Y, size = mismatch_size, replace = TRUE)
 }
 return(data_obs)
}

################################################################################
# Estimation and Variance Framework (Section 3)
################################################################################
# Estimator 1 & 2: Plain Adjustment EM with Analytic SE (Louis' Method)
fit_plain_em <- function(tab, alpha, tol = 1e-6, max_iter = 1000) {
 n_total <- sum(tab)
 n_cells <- length(tab)

 p_sat <- tab / n_total
 p_ind <- outer(rowSums(p_sat), colSums(p_sat))
 eps <- .Machine$double.eps

 # EM Algorithm Loop
 for (iter in 1:max_iter) {
  numer <- p_sat * (1 - alpha)
  denom <- numer + (p_ind * alpha)
  weights <- numer / (denom + eps)

  counts_sat <- tab * weights
  if (sum(counts_sat) > 0) p_sat <- counts_sat / sum(counts_sat)

  counts_ind <- tab * (1 - weights)
  sum_ind <- sum(counts_ind)
  if (sum_ind > 1e-12) p_ind <- outer(rowSums(counts_ind)/sum_ind, colSums(counts_ind)/sum_ind)

  if (iter > 1 && max(abs(p_sat - p_sat_prev)) < tol) break
  p_sat_prev <- p_sat
 }

 # Asymptotic Variance Extraction (Louis 1982 via Slawski et al., 2025)
 counts_vec <- c(t(tab))
 p_sat_vec  <- c(t(p_sat))
 p_ind_vec  <- c(t(p_ind))

 denom_sq <- ((1 - alpha) * p_sat_vec + alpha * p_ind_vec)^2
 hess0 <- counts_vec * (1 - alpha)^2 / (denom_sq + eps)

 ones <- rep(1, n_cells)
 term_outer <- outer(ones / n_cells, hess0) + outer(hess0, ones / n_cells)

 hess <- -term_outer + mean(hess0) / n_cells
 diag(hess) <- diag(hess) + hess0

 svd_hess <- svd(hess)
 rank_idx <- seq_len(max(1, n_cells - 1))
 inv_d <- 1 / sqrt(svd_hess$d[rank_idx])
 root_inv <- scale(svd_hess$u[, rank_idx, drop = FALSE], center = FALSE, scale = 1/inv_d)
 vcov_phat <- tcrossprod(root_inv)

 se_mat <- matrix(sqrt(diag(vcov_phat)), nrow = nrow(tab), byrow = TRUE)

 return(list(theta = p_sat, se = se_mat))
}

# Estimator 3: Proposed Three-Component EM (Section 3.2.2)
fit_proposed_em <- function(C, data_obs, s, alpha, tol = 1e-6, max_iter = 1000) {
 # Calculate Fixed Mixing Weights
 w1 <- (1 - alpha) * (1 - s)
 w2 <- (1 - alpha) * s
 w3 <- alpha
 N <- sum(C)
 eps <- .Machine$double.eps

 # Pre-compute fixed noise margins from observed data
 P_ind <- outer(prop.table(table(data_obs$X)), prop.table(table(data_obs$Y)))

 P_CI <- matrix(0, nrow = 2, ncol = 2)
 for (z in unique(data_obs$Z)) {
  sub_data <- subset(data_obs, Z == z)
  p_z <- nrow(sub_data) / N
  if (nrow(sub_data) > 0) {
   P_CI <- P_CI + outer(prop.table(table(factor(sub_data$X, levels=c(0,1)))),
                        prop.table(table(factor(sub_data$Y, levels=c(0,1))))) * p_z
  }
 }

 theta <- prop.table(C)

 # EM Algorithm Loop
 for (iter in 1:max_iter) {
  # E-Step: Posterior allocation weights
  denom <- (w1 * theta) + (w2 * P_CI) + (w3 * P_ind) + eps
  tau <- (w1 * theta) / denom

  # M-Step: Update parameters
  C_tilde <- C * tau
  theta_new <- C_tilde / sum(C_tilde)

  if (max(abs(theta_new - theta)) < tol) {
   theta <- theta_new
   break
  }
  theta <- theta_new
 }
 return(list(theta = theta))
}

# Variance Estimation: Non-Parametric Bootstrap (Section 3.2.3.2)
bootstrap_proposed_nonparametric <- function(data_obs, s, alpha, B = 200) {
 N <- nrow(data_obs)
 boot_ests <- matrix(NA, nrow = B, ncol = 4)

 for (b in 1:B) {
  # 1. Resample directly from the observed data (Preserves conditional Z structure)
  boot_idx <- sample(1:N, size = N, replace = TRUE)
  boot_obs <- data_obs[boot_idx, ]

  # 2. Extract bootstrapped counts
  C_boot <- table(factor(boot_obs$X, levels=c(0,1)), factor(boot_obs$Y, levels=c(0,1)))

  # 3. Re-run EM (The function internally recomputes P_ind and P_CI on boot_obs)
  boot_fit <- fit_proposed_em(C_boot, boot_obs, s, alpha)

  # Store flattened estimates
  boot_ests[b, ] <- as.vector(t(boot_fit$theta))
 }

 # Calculate empirical standard errors
 se_vec <- apply(boot_ests, 2, sd)
 return(matrix(se_vec, nrow = 2, byrow = TRUE))
}

################################################################################
# Core Simulation Loop (Section 4.3)
################################################################################
run_sim_iteration <- function(iter, s, alpha, N = 5000, B = 200) {

 # A. Analytical DGP Target
 true_11 <- (0.3 * 0.7 * (0.8 + 0.2 * 0.7)) +
  (0.4 * 0.5 * (0.9 + 0.1 * 0.5)) +
  (0.3 * 0.3 * (0.7 + 0.3 * 0.3))

 # B. Latent Truth & Degradation
 data_true <- generate_population(N)
 P_oracle <- prop.table(table(factor(data_true$X, levels=c(0,1)), factor(data_true$Y, levels=c(0,1))))

 data_obs <- perturb_data(data_true, s, alpha)
 C_obs <- table(factor(data_obs$X, levels=c(0,1)), factor(data_obs$Y, levels=c(0,1)))
 P_naive <- prop.table(C_obs)

 # C. Estimation, Inference, and Time Tracking

 # 0. Oracle Estimator
 t0 <- Sys.time()
 Oracle_11 <- P_oracle[2,2]
 Oracle_SE <- sqrt(Oracle_11 * (1 - Oracle_11) / N)
 time_oracle <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 # 1. Naive Estimator
 t0 <- Sys.time()
 Naive_11 <- P_naive[2,2]
 Naive_SE <- sqrt(Naive_11 * (1 - Naive_11) / N)
 time_naive <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 # 2. Plain Adjustment
 t0 <- Sys.time()
 fit_plain <- fit_plain_em(C_obs, alpha)
 time_plain <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 # 3. Plain Adjustment Gamma
 t0 <- Sys.time()
 gamma <- alpha + s - (alpha * s)
 fit_gamma <- fit_plain_em(C_obs, gamma)
 time_gamma <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 # 4. Proposed Three-Component Framework
 t0 <- Sys.time()
 fit_prop <- fit_proposed_em(C_obs, data_obs, s, alpha)
 se_prop <- bootstrap_proposed_nonparametric(data_obs, s, alpha, B = B)
 time_prop <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

 # Helper to compute 95% Wald interval coverage
 calc_coverage <- function(est, se, truth) {
  ifelse(truth >= (est - 1.96 * se) & truth <= (est + 1.96 * se), 1, 0)
 }

 res <- data.frame(
  iter = iter, s = s, alpha = alpha, True_11 = true_11,

  Oracle_11 = Oracle_11, Oracle_SE = Oracle_SE, Oracle_Cov = calc_coverage(Oracle_11, Oracle_SE, true_11),
  Naive_11 = Naive_11, Naive_SE = Naive_SE, Naive_Cov = calc_coverage(Naive_11, Naive_SE, true_11),
  Plain_11 = fit_plain$theta[2,2], Plain_SE = fit_plain$se[2,2], Plain_Cov = calc_coverage(fit_plain$theta[2,2], fit_plain$se[2,2], true_11),
  Gamma_11 = fit_gamma$theta[2,2], Gamma_SE = fit_gamma$se[2,2], Gamma_Cov = calc_coverage(fit_gamma$theta[2,2], fit_gamma$se[2,2], true_11),
  Prop_11 = fit_prop$theta[2,2], Prop_SE = se_prop[2,2], Prop_Cov = calc_coverage(fit_prop$theta[2,2], se_prop[2,2], true_11),

  # Run Times
  Time_Oracle = time_oracle,
  Time_Naive = time_naive,
  Time_Plain = time_plain,
  Time_Gamma = time_gamma,
  Time_Prop = time_prop
 )

 return(res)
}

################################################################################
# Windows Parallel Execution Grid (Section 4.1)
################################################################################
# Perturbation parameters grid
s_grid <- c(0.05, 0.15, 0.25)
alpha_grid <- c(0.05, 0.10, 0.20)
params <- expand.grid(s = s_grid, alpha = alpha_grid)

# Iteration counts
M <- 1000 #1000
B <- 200  #200

# Initialize Socket Cluster for Windows
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterEvalQ(cl, { library(dplyr); library(tidyr) })

# Export all required objects to the Windows worker environments
clusterExport(cl, c("run_sim_iteration", "generate_population", "perturb_data",
                    "fit_plain_em", "fit_proposed_em",
                    "bootstrap_proposed_nonparametric"))

results_list <- list()

# Starting Full Simulation Grid
for (i in 1:nrow(params)) {
 s_current <- params$s[i]
 alpha_current <- params$alpha[i]

 cat(sprintf("Running grid point %d/%d: s = %.2f, alpha = %.2f\n",
             i, nrow(params), s_current, alpha_current))

 # Pass s_current, alpha_current, and B into the parLapply environment
 sim_out <- parLapply(cl, 1:M, function(iter, s_val, alpha_val, B_val) {
  run_sim_iteration(iter, s = s_val, alpha = alpha_val, N = 5000, B = B_val)
 }, s_val = s_current, alpha_val = alpha_current, B_val = B)

 results_list[[i]] <- bind_rows(sim_out)
}

# Close workers
stopCluster(cl)

################################################################################
# Summarize Results: Bias, RMSE, and Coverage
################################################################################
all_results <- bind_rows(results_list)

summary_table <- all_results %>%
 group_by(s, alpha) %>%
 summarize(
  Truth_Avg = mean(True_11),

  Bias_Oracle= mean(Oracle_11 - True_11), Bias_Naive = mean(Naive_11 - True_11),
  Bias_Plain = mean(Plain_11 - True_11), Bias_Gamma = mean(Gamma_11 - True_11),
  Bias_Prop = mean(Prop_11 - True_11),

  RMSE_Oracle= sqrt(mean((Oracle_11 - True_11)^2)),
  RMSE_Naive = sqrt(mean((Naive_11 - True_11)^2)),
  RMSE_Plain = sqrt(mean((Plain_11 - True_11)^2)),
  RMSE_Gamma = sqrt(mean((Gamma_11 - True_11)^2)),
  RMSE_Prop = sqrt(mean((Prop_11 - True_11)^2)),

  Cov_Oracle = mean(Oracle_Cov), Cov_Naive = mean(Naive_Cov),
  Cov_Plain = mean(Plain_Cov),
  Cov_Gamma  = mean(Gamma_Cov), Cov_Prop = mean(Prop_Cov),

  # Average Run Time per Iteration (Seconds)
  Time_Oracle = mean(Time_Oracle),
  Time_Naive = mean(Time_Naive),
  Time_Plain = mean(Time_Plain),
  Time_Gamma = mean(Time_Gamma),
  Time_Prop  = mean(Time_Prop),

  .groups = 'drop'
 )

print(round(summary_table, 4))

################################################################################
# Generate Plots for Simulation Results
################################################################################
plot_data <- summary_table %>%
 select(-Truth_Avg) %>%
 pivot_longer(
  cols = starts_with("Bias_") | starts_with("RMSE_") | starts_with("Cov_"),
  names_to = c("Metric", "Estimator"),
  names_sep = "_"
 ) %>%
 mutate(
  Metric = factor(Metric,
                  levels = c("Bias", "RMSE", "Cov"),
                  labels = c("Absolute Bias", "RMSE", "Coverage Probability")),
  Estimator = factor(Estimator,
                     levels = c("Oracle", "Naive", "Plain", "Gamma", "Prop"),
                     labels = c("Oracle", "Naive", "Plain Adj.", "Gamma Adj.",
                                "Proposed")),
  alpha_factor = factor(sprintf("%.2f", alpha))
 )

ref_lines <- data.frame(
 Metric = factor(c("Absolute Bias", "RMSE", "Coverage Probability"),
                 levels = c("Absolute Bias", "RMSE", "Coverage Probability")),
 yintercept = c(0, 0, 0.95)
)

pd <- position_dodge(width = 0.05)

my_colors <- c("Oracle" = "grey50", "Naive" = "#D55E00",
               "Plain Adj." = "#0072B2", "Gamma Adj." = "#009E73",
               "Proposed" = "#CC79A7")
my_shapes <- c("Oracle" = 15, "Naive" = 16, "Plain Adj." = 17,
               "Gamma Adj." = 18, "Proposed" = 19)
my_lines  <- c("Oracle" = "solid", "Naive" = "solid",
               "Plain Adj." = "dashed", "Gamma Adj." = "dotdash",
               "Proposed" = "solid")
my_widths <- c("Oracle" = 3, "Naive" = 0.8, "Plain Adj." = 0.8,
               "Gamma Adj." = 0.8, "Proposed" = 0.8)
my_alphas <- c("Oracle" = 0.3, "Naive" = 1, "Plain Adj." = 1,
               "Gamma Adj." = 1, "Proposed" = 1)

p <- ggplot(plot_data, aes(x = alpha_factor, y = value,
                           color = Estimator, shape = Estimator,
                           linetype = Estimator, linewidth = Estimator,
                           alpha = Estimator,
                           group = Estimator)) +

 geom_hline(data = ref_lines, aes(yintercept = yintercept),
            linetype = "solid", color = "black", alpha = 0.3, linewidth = 0.8) +

 geom_line(position = pd) +
 geom_point(position = pd, size = 3.5, stroke = 1) +

 facet_grid(Metric ~ s, scales = "free_y",
            labeller = labeller(s = function(x) paste("Swap Rate (s) =", x))) +

 labs(x = expression(paste("Mismatch Rate (", alpha, ")")), y = NULL) +
 scale_color_manual(values = my_colors) +
 scale_shape_manual(values = my_shapes) +
 scale_linetype_manual(values = my_lines) +
 scale_linewidth_manual(values = my_widths) +
 scale_alpha_manual(values = my_alphas) +

 scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +

 theme_bw(base_size = 14) +
 theme(
  strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
  strip.text = element_text(face = "bold", size = 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", linewidth = 1),
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.key.width = unit(3, "cm"),
  legend.text = element_text(size = 12),
  axis.title.x = element_text(margin = margin(t = 15), face = "bold")
 )

print(p)
ggsave("Simulation_Results.pdf", plot = p, width = 11, height = 8.5, dpi = 600)

################################################################################
# Generate Table for Simulation Results
################################################################################
formatted_data <- summary_table %>%
 select(s, alpha,
        Bias_Oracle, Bias_Naive, Bias_Plain, Bias_Gamma, Bias_Prop,
        RMSE_Oracle, RMSE_Naive, RMSE_Plain, RMSE_Gamma, RMSE_Prop,
        Cov_Oracle,  Cov_Naive,  Cov_Plain,  Cov_Gamma,  Cov_Prop,
        Time_Oracle, Time_Naive, Time_Plain, Time_Gamma, Time_Prop) %>%
 mutate(across(starts_with("Bias"), ~ sprintf("%.3f", .)),
        across(starts_with("RMSE"), ~ sprintf("%.3f", .)),
        across(starts_with("Cov"),  ~ sprintf("%.2f", .)),
        across(starts_with("Time"), ~ sprintf("%.4f", .))) %>%
 arrange(s, alpha)

latex_str <- c(
 "\\begin{table}[ht]",
 "\\centering",
 "\\caption{Simulation results: Statistical performance and average computational run times for the target estimand $P(X=1, Y=1)$ across all estimators.}",
 "\\label{tab:sim_results_full_with_time}",
 "% Resize to fit the page margins automatically",
 "\\resizebox{\\textwidth}{!}{",
 "\\begin{tabular}{cc|ccccc|ccccc|ccccc}",
 "\\toprule",
 "\\multicolumn{17}{c}{\\textbf{Panel A: Statistical Inference Performance}} \\\\",
 "\\midrule",
 "\\multicolumn{2}{c}{Design} & \\multicolumn{5}{c}{Absolute Bias} & \\multicolumn{5}{c}{RMSE} & \\multicolumn{5}{c}{Coverage (95\\%)} \\\\",
 "\\cmidrule(lr){1-2} \\cmidrule(lr){3-7} \\cmidrule(lr){8-12} \\cmidrule(lr){13-17}",
 "$s$ & $\\alpha$ & Oracle & Naive & Plain & Gamma & Prop & Oracle & Naive & Plain & Gamma & Prop & Oracle & Naive & Plain & Gamma & Prop \\\\",
 "\\midrule"
)

current_s <- ""
for (i in 1:nrow(formatted_data)) {
 row <- formatted_data[i, ]

 if (row$s != current_s) {
  if (i > 1) latex_str <- c(latex_str, "\\midrule")
  s_col <- sprintf("\\multirow{3}{*}{%.2f}", as.numeric(row$s))
  current_s <- row$s
 } else {
  s_col <- ""
 }

 row_str <- sprintf("%s & %.2f & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                    s_col, as.numeric(row$alpha),
                    row$Bias_Oracle, row$Bias_Naive, row$Bias_Plain, row$Bias_Gamma, row$Bias_Prop,
                    row$RMSE_Oracle, row$RMSE_Naive, row$RMSE_Plain, row$RMSE_Gamma, row$RMSE_Prop,
                    row$Cov_Oracle, row$Cov_Naive, row$Cov_Plain, row$Cov_Gamma, row$Cov_Prop)

 latex_str <- c(latex_str, row_str)
}

latex_str <- c(
 latex_str,
 "\\midrule",
 "\\midrule",
 "\\multicolumn{17}{c}{\\textbf{Panel B: Average Computational Run Time (Seconds per Iteration)}} \\\\",
 "\\midrule",
 "$s$ & $\\alpha$ & \\multicolumn{3}{c}{Oracle} & \\multicolumn{3}{c}{Naive} & \\multicolumn{3}{c}{Plain Adj.} & \\multicolumn{3}{c}{Gamma Adj.} & \\multicolumn{3}{c}{Proposed} \\\\",
 "\\midrule"
)

current_s <- ""
for (i in 1:nrow(formatted_data)) {
 row <- formatted_data[i, ]

 if (row$s != current_s) {
  if (i > 1) latex_str <- c(latex_str, "\\midrule")
  s_col <- sprintf("\\multirow{3}{*}{%.2f}", as.numeric(row$s))
  current_s <- row$s
 } else {
  s_col <- ""
 }

 row_str <- sprintf("%s & %.2f & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} \\\\",
                    s_col, as.numeric(row$alpha),
                    row$Time_Oracle, row$Time_Naive, row$Time_Plain, row$Time_Gamma, row$Time_Prop)

 latex_str <- c(latex_str, row_str)
}

latex_str <- c(
 latex_str,
 "\\bottomrule",
 "\\end{tabular}",
 "}",
 "\\end{table}"
)

cat(paste(latex_str, collapse = "\n"))
