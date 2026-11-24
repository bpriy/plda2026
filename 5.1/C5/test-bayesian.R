### data generation
ls()

mismatches_fun <- function(gammaval){
 yperm <- y
 
 #' generate mismatches (w/ z = X)
 m <- 1 - rbinom(n, prob = plogis(X %*% gammaval), size = 1) 
 mmr <- sum(m)/n
 
 yperm <- y
 # only shuffle the mismatched data
 shuffled_ix <- sample(which(m == 1))
 yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], 
                               shuffled_ix[1])]
 
 return(list(yperm = yperm, mmr = mmr, m = m))
}

#' Set random seed 
set.seed(124)

#' Generate correctly linked data (X)
n <- 1000 # number of observations
p <- 1 # number of covariates
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X <- cbind(1, X) # intercept + covariates

true_beta <- c(1, -1) # true beta
true_gamma <- c(2.5, 4.5) # true gamma
true_sigma <- 0.25 # true sigma

y <- rnorm(n, X %*% true_beta, true_sigma) # true y

#' Generate mismatched y
mismatches <- mismatches_fun(true_gamma)
yperm <- mismatches$yperm

Delta <- X 

#  use yperm instead of y to extract eta 
eta <- quantile(yperm, seq(from = 0.02, to = 0.98, length = 50))

### STAN
library(rstan)

#  Removed the double commas from the lists
init_values <- list(
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma, tau = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma, tau = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma, tau = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma, tau = true_sigma)
)

#  Change to x = X[,2] and z = X[,2] so the covariate is passed, not the intercept
fit <- stan(file = 'ext1-b.stan', 
            data = list(n = 1000, m = 50, x = X[,2], y = yperm, z = X[,2], eta = as.numeric(eta)), 
            pars = c("beta0", "beta1", "gamma0", "gamma1", "sigma", "tau", "alpha"), 
            init = init_values,
            seed = 12345) # seed for reproducibility

summary_fit <- summary(fit)

#  Print parameters to check the results
print(fit, pars = c("beta0", "beta1", "gamma0", "gamma1", "sigma", "tau"))

### JAGS
library(rjags)

#  Changed X[,1] to X[,2] so we pass the covariate instead of the intercept column.
#  Added alpha_prior to avoid using rep(1, m) syntax in JAGS.
#  Removed sigmainv, tauinv, gamma0, and gamma1 from the dat list since they are now estimated.
dat <- list(x = X[,2], z = X[,2], y = yperm, n = n, m = length(eta), 
            eta = c(0, eta), alpha_prior = rep(1, length(eta)))

jags_inits <- list(beta0 = true_beta[1],
                   beta1 = true_beta[2],
                   .RNG.name = "base::Wichmann-Hill", # RNG
                   .RNG.seed = 12345)                 # seed

jm <- jags.model(file = "ext1-jags.txt", data = dat,
                 n.chains = 1, n.adapt=1000, quiet=FALSE,
                 inits = jags_inits)

#  Added to run the sampling and display results
update(jm, 1000)
#  Added gamma0, gamma1, sigma, and tau to the tracked variables
samples <- coda.samples(jm, 
                        variable.names = c("beta0", "beta1", "gamma0", "gamma1",
                                           "sigma", "tau"), 
                        n.iter = 2000)
summary(samples)

################################################################################
## POSTERIOR SUMMARIES
################################################################################
# parameters of interest for Formulation 2
params <- c("beta0", "beta1", "gamma0", "gamma1", "sigma", "tau")

jags_sum <- summary(samples)

# Combine Mean, SD, and the quantiles
jags_table <- cbind(
 mean = jags_sum$statistics[params, "Mean"],
 sd   = jags_sum$statistics[params, "SD"],
 jags_sum$quantiles[params, c("2.5%", "25%", "50%", "75%", "97.5%")]
)

print(round(jags_table, 3))

# Extract summary specifying the exact probabilities needed
stan_sum <- summary(fit, pars = params, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
# Subset to the columns that match the supplement format
stan_table <- stan_sum[, c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")]
print(round(stan_table, 3))

################################################################################
## TRACEPLOTS
################################################################################
# Export JAGS traceplots
pdf("jags_traceplot_2.pdf", width = 7, height = 7)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1)) 
for (p in params) {
 chain_matrix <- as.matrix(samples[, p, drop = FALSE])
 matplot(chain_matrix, 
         type = "l", 
         lty = 1,
         col = "black", 
         main = p, 
         xlab = "Iterations", 
         ylab = "")
}
par(mfrow = c(1, 1)) # Reset graphic parameters
dev.off()

# Export STAN traceplots
pdf("stan_traceplot_2.pdf", width = 7, height = 7)
stan_traceplot <- rstan::traceplot(fit, pars = params, inc_warmup = FALSE)
print(stan_traceplot)
dev.off()