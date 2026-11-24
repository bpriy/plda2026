library(tictoc)
source("pldamixture_glm.R")

################################################################################
## PMM FUNCTION
################################################################################
project_MNAR_PMM <- function(S1_params, S0_moments, phi) {
 # Extract S=1 (Adjusted) Parameters
 b_Y1    <- S1_params$slope_Y
 b0_Y1   <- S1_params$int_Y
 b_X1    <- S1_params$slope_X
 b0_X1   <- S1_params$int_X
 var_Y1  <- S1_params$var_Y
 var_X1  <- S1_params$var_X
 rho_XY1 <- S1_params$rho_XY
 
 # Extract S=0 Summary Moments
 b_X0  <- S0_moments$slope_X
 b0_X0 <- S0_moments$int_X
 
 # Calculate the sensitivity scaling multiplier g(phi)
 g_phi <- (phi + (1 - phi) * rho_XY1) / (phi * rho_XY1 + (1 - phi))
 
 # Recover Unobserved Target Slope
 b_Y0  <- b_Y1 + g_phi * sqrt(var_Y1 / var_X1) * (b_X0 - b_X1)
 
 # Recover Unobserved Target Intercept (West et al., 2021 Eq 11)
 b0_Y0 <- b0_Y1 + g_phi * sqrt(var_Y1 / var_X1) * (b0_X0 - b0_X1)
 
 return(list(
  recovered_slope = b_Y0,
  recovered_int   = b0_Y0,
  g_multiplier    = g_phi
 ))
}

################################################################################
## GENERATE LINKED DATA
################################################################################
# Initialize Profiling
tic.clearlog()
tic("Total Simulation Time", log = TRUE)

tic("Data Generating Process (MNAR) & Introduce Linkage Error")
set.seed(999)
n <- 10000

# Generate Primary Covariate Z
Z <- rnorm(n, 0, 1)

# Generate Auxiliary Covariates
A1 <- rnorm(n, 0.5 * Z, 1)
A2 <- rnorm(n, -0.5 * Z, 1)

# Generate Latent Outcome Y
Y <- -1 + 2 * Z + 1.5 * A1 - 1.2 * A2 + rnorm(n, 0, 1)

# Modified Non-Ignorable Selection Mechanism
# By balancing the main effects (0.5*Z and -0.2*Y), we minimize the between-group 
# differences in means, making the overall slope bounded between S=1 and S=0, 
# while keeping the Z*Y interaction (0.8) to ensure the slopes differ.
delta <- rbinom(n, size = 1, prob = plogis(0.5 + 0.5*Z - 0.2*Y + 0.8*Z*Y))

Ycomp <- Y[delta == 1]; Zcomp <- Z[delta == 1]
A1comp <- A1[delta == 1]; A2comp <- A2[delta == 1]

Ymis <- Y[delta == 0]; Zmis <- Z[delta == 0]
A1mis <- A1[delta == 0]; A2mis <- A2[delta == 0]

# Introduce 15% Mismatch Error into Linked Covariates
m.rate <- 0.15
n_comp <- length(Zcomp)
idx_mismatch <- sample(1:n_comp, size = floor(n_comp * m.rate))
shuffled_idx <- sample(idx_mismatch)

Zcomp_linked <- Zcomp
A1comp_linked <- A1comp
A2comp_linked <- A2comp

Zcomp_linked[idx_mismatch] <- Zcomp[shuffled_idx]
A1comp_linked[idx_mismatch] <- A1comp[shuffled_idx]
A2comp_linked[idx_mismatch] <- A2comp[shuffled_idx]

M_Z <- matrix(rnorm(n_comp, mean = ifelse(1:n_comp %in% idx_mismatch, -1, 1)), 
              ncol = 1)
toc(log = TRUE) # End Data Generation

################################################################################
## SIMULATION STUDY
################################################################################
tic("Step 1: Proxy Construction & Adjustment (S=1)")
tic("Step 1a: Auxiliary Weight Estimation via Mixture Model")
fit_YZA <- glmMixture(
 x = cbind(1, Zcomp_linked, A1comp_linked, A2comp_linked),
 y = Ycomp, family = "gaussian", z = cbind(1, M_Z)
)
beta_a <- fit_YZA$coefficients[c("A1comp_linked", "A2comp_linked")]
toc(log = TRUE)

# Construct proxy X
Xcomp_linked <- A1comp_linked * beta_a[1] + A2comp_linked * beta_a[2]
Xmis <- A1mis * beta_a[1] + A2mis * beta_a[2]

tic("Step 1b: Parameter Estimation", log = TRUE)
fit_Y_Z <- glmMixture(x = cbind(1, Zcomp_linked), y = Ycomp, 
                      family = "gaussian", z = cbind(1, M_Z))
fit_X_Z <- glmMixture(x = cbind(1, Zcomp_linked), y = Xcomp_linked, 
                      family = "gaussian", z = cbind(1, M_Z))
fit_X_YZ <- glmMixture(x = cbind(1, Ycomp, Zcomp_linked), y = Xcomp_linked, 
                       family = "gaussian", z = cbind(1, M_Z))
toc(log = TRUE)

# Extract variances to compute conditional correlation (rho)
var_Y1 <- fit_Y_Z$dispersion
var_X1 <- fit_X_Z$dispersion
cov_XY1 <- fit_X_YZ$coefficients["Ycomp"] * var_Y1
rho_XY1 <- cov_XY1 / sqrt(var_X1 * var_Y1)

# Store S=1 parameters
adj_S1_params <- list(
 slope_Y = fit_Y_Z$coefficients["Zcomp_linked"],
 int_Y   = fit_Y_Z$coefficients[[1]], 
 slope_X = fit_X_Z$coefficients["Zcomp_linked"],
 int_X   = fit_X_Z$coefficients[[1]],
 var_Y   = var_Y1,
 var_X   = var_X1,
 rho_XY  = rho_XY1
)

toc(log = TRUE) # End Step 1


tic("Steps 2 & 3: Summary Moments (S=0) & MNAR PMM Projection")
# Extract S=0 Summary Moments 
lm_X_mis <- lm(Xmis ~ Zmis)
S0_moments <- list(
 slope_X = coef(lm_X_mis)["Zmis"],
 int_X   = coef(lm_X_mis)["(Intercept)"]
)

phi_grid <- c(0, 0.5, 1)
results_proj <- list()
for (phi in phi_grid) {
 results_proj[[as.character(phi)]] <- project_MNAR_PMM(
  S1_params = adj_S1_params,
  S0_moments = S0_moments,
  phi = phi
 )
}

toc(log = TRUE) # End Steps 2 & 3


tic("Benchmarking & Performance Evaluation")
# Oracle (Ground Truth S=0)
oracle_beta_Y_Z_1 <- coef(lm(Ycomp ~ Zcomp))["Zcomp"]
oracle_beta_Y_Z_0 <- coef(lm(Ymis ~ Zmis))["Zmis"]

# Doubly Naive (Ignores mismatches and ignores selection bias)
lm_doubly_naive <- lm(Ycomp ~ Zcomp_linked)
doubly_naive_slope <- coef(lm_doubly_naive)["Zcomp_linked"]

# Naive Selection-Adjusted (Ignores mismatches, accounts for selection bias)
# Re-estimate proxy weights using naive OLS (ignoring mismatch)
naive_fit_YZA <- lm(Ycomp ~ Zcomp_linked + A1comp_linked + A2comp_linked)
naive_beta_a <- coef(naive_fit_YZA)[c("A1comp_linked", "A2comp_linked")]

# Construct naive proxy
naive_Xcomp_linked <- A1comp_linked * naive_beta_a[1] + A2comp_linked * 
 naive_beta_a[2]
naive_Xmis <- A1mis * naive_beta_a[1] + A2mis * naive_beta_a[2]

# Get naive S=1 parameters
naive_fit_Y_Z <- lm(Ycomp ~ Zcomp_linked)
naive_fit_X_Z <- lm(naive_Xcomp_linked ~ Zcomp_linked)
naive_fit_X_YZ <- lm(naive_Xcomp_linked ~ Ycomp + Zcomp_linked)

naive_var_Y1 <- sigma(naive_fit_Y_Z)^2
naive_var_X1 <- sigma(naive_fit_X_Z)^2
naive_cov_XY1 <- coef(naive_fit_X_YZ)["Ycomp"] * naive_var_Y1
naive_rho_XY1 <- naive_cov_XY1 / sqrt(naive_var_X1 * naive_var_Y1)

naive_S1_params <- list(
 slope_Y = coef(naive_fit_Y_Z)["Zcomp_linked"],
 int_Y   = coef(naive_fit_Y_Z)["(Intercept)"],
 slope_X = coef(naive_fit_X_Z)["Zcomp_linked"],
 int_X   = coef(naive_fit_X_Z)["(Intercept)"],
 var_Y   = naive_var_Y1,
 var_X   = naive_var_X1,
 rho_XY  = naive_rho_XY1
)

lm_naive_X_mis <- lm(naive_Xmis ~ Zmis)
naive_S0_moments <- list(
 slope_X = coef(lm_naive_X_mis)["Zmis"],
 int_X   = coef(lm_naive_X_mis)["(Intercept)"]
)

# Test naive selection-adjustment at phi = 0.5 baseline
naive_sel_adj_proj <- project_MNAR_PMM(S1_params = naive_S1_params, 
                                       S0_moments = naive_S0_moments, 
                                       phi = 0.5)

# Naive Mismatch-Adjusted (Accounts for mismatches via glmMixture, 
# ignores selection bias)
naive_mis_adj_slope <- adj_S1_params$slope_Y

# Get values for Overall Slope 
# Get the empirical Oracle overall slope (Full Population Ground Truth)
oracle_overall <- coef(lm(Y ~ Z))["Z"] 

# Calculate the marginal non-selection and selection rates
pr_delta_0 <- mean(delta == 0)
pr_delta_1 <- mean(delta == 1)

# Compute Overall Slopes 
# Helper function to compute Overall Slope via Law of Total Covariance
compute_LTC_slope <- function(b1, b0, int1, int0) {
 # Conditional Covariances
 cov_YZ1 <- b1 * var(Zcomp_linked)
 cov_YZ0 <- b0 * var(Zmis)
 
 # Conditional Means of Y (using linear equation E[Y|S] = int + b * E[Z|S])
 mu_Y1 <- int1 + b1 * mean(Zcomp_linked)
 mu_Y0 <- int0 + b0 * mean(Zmis)
 
 # Law of Total Covariance
 total_cov <- (pr_delta_1 * cov_YZ1) + (pr_delta_0 * cov_YZ0) + 
  (pr_delta_1 * pr_delta_0 * (mu_Y1 - mu_Y0) * (mean(Zcomp_linked) - mean(Zmis)))
 
 # Law of Total Variance
 total_var <- (pr_delta_1 * var(Zcomp_linked)) + (pr_delta_0 * var(Zmis)) + 
  (pr_delta_1 * pr_delta_0 * (mean(Zcomp_linked) - mean(Zmis))^2)
 
 return(total_cov / total_var)
}

# Doubly Naive Overall (Assumes S=0 is identical to S=1)
b_dn <- coef(lm_doubly_naive)["Zcomp_linked"]
int_dn <- coef(lm_doubly_naive)["(Intercept)"]
overall_doubly_naive <- compute_LTC_slope(b_dn, b_dn, int_dn, int_dn)

# Naive Mismatch-Adjusted Only (Assumes S=0 is identical to S=1)
overall_naive_mis_adj <- compute_LTC_slope(adj_S1_params$slope_Y, 
                                           adj_S1_params$slope_Y, 
                                           adj_S1_params$int_Y,
                                           adj_S1_params$int_Y)

# Naive Selection-Adjusted Only (Mix naive S=1 with naive-projected S=0)
overall_naive_sel_adj <- compute_LTC_slope(naive_S1_params$slope_Y, 
                                           naive_sel_adj_proj$recovered_slope, 
                                           naive_S1_params$int_Y, 
                                           naive_sel_adj_proj$recovered_int)

# Proposed Adjusted Methods
overall_phi_0   <- compute_LTC_slope(adj_S1_params$slope_Y,
                                     results_proj[["0"]]$recovered_slope, 
                                     adj_S1_params$int_Y, 
                                     results_proj[["0"]]$recovered_int)
overall_phi_0.5 <- compute_LTC_slope(adj_S1_params$slope_Y, 
                                     results_proj[["0.5"]]$recovered_slope, 
                                     adj_S1_params$int_Y, 
                                     results_proj[["0.5"]]$recovered_int)
overall_phi_1   <- compute_LTC_slope(adj_S1_params$slope_Y, 
                                     results_proj[["1"]]$recovered_slope, 
                                     adj_S1_params$int_Y, 
                                     results_proj[["1"]]$recovered_int)

# Compile Results
results_df <- data.frame(
 Parameter = c("Slope (S=1)", "Slope (S=0)", "Overall Slope"),
 Oracle = c(oracle_beta_Y_Z_1, oracle_beta_Y_Z_0, oracle_overall),
 Doubly_Naive = c(doubly_naive_slope, doubly_naive_slope, overall_doubly_naive),
 Naive_Sel_Adj_Only = c(naive_S1_params$slope_Y, 
                     naive_sel_adj_proj$recovered_slope, overall_naive_sel_adj),
 Naive_Mis_Adj_Only = c(naive_mis_adj_slope, 
                        naive_mis_adj_slope, overall_naive_mis_adj),
 Proposed_phi_0 = c(adj_S1_params$slope_Y, 
                    results_proj[["0"]]$recovered_slope, overall_phi_0),
 Proposed_phi_0.5 = c(adj_S1_params$slope_Y, 
                      results_proj[["0.5"]]$recovered_slope, overall_phi_0.5),
 Proposed_phi_1 = c(adj_S1_params$slope_Y, 
                    results_proj[["1"]]$recovered_slope, overall_phi_1)
)

knitr::kable(results_df, digits = 3, row.names = FALSE,
             col.names = c("Parameter", "Oracle", "Doubly Naïve", 
                           "Naïve (Sel. Adj. Only)", "Naïve (Mis. Adj. Only)",
                           "Proposed (phi=0)", "(phi=0.5)", "(phi=1)"))

toc(log = TRUE) # End Benchmarking

toc(log = TRUE) # End Total Simulation Time

# Print Profile
unlist(tic.log())