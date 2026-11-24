### data generation
rm(list=ls())

mismatches_fun <- function(gammaval){
 yperm <- y
 
 #' generate mismatches (w/ z = X)
 m <- 1 - rbinom(n, prob = plogis(X %*% gammaval), size = 1) 
 mmr <- sum(m)/n
 
 yperm <- y
 # only shuffle the mismatched data
 shuffled_ix <- sample(which(m == 1))
 if(length(shuffled_ix) > 1) {
  yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], 
                                shuffled_ix[1])]
 }
 
 return(list(yperm = yperm, mmr = mmr, m = m))
}

#' Set random seed 
set.seed(124)

#' Generate correctly linked data (X)
n <- 150 # reduced n to 150 for tractability (originally 1000)
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

### STAN
library(rstan)

# Prepare initial values list for Stan
init_values <- list(
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma),
 list(gamma0 = true_gamma[1], gamma1 = true_gamma[2], beta0 = true_beta[1], beta1 = true_beta[2], sigma = true_sigma)
)

# Compile and run the Stan model
# 'eta' and 'm' are removed from the data list
fit <- stan(file = 'ext1o-b.stan', 
            data = list(n = n, x = X[,2], y = yperm, z = X[,2]), 
            pars = c("beta0", "beta1", "gamma0", "gamma1", "sigma"), 
            init = init_values,
            iter = 1000,     # reduce iter from 2k
            warmup = 500,    # set warmup to half
            chains = 4,
            cores = 4) # Run chains in parallel for speed

# Display results
print(fit)
traceplot(fit)