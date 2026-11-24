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
 if (length(shuffled_ix) > 1) {
  yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], 
                                shuffled_ix[1])]
 }
 
 return(list(yperm = yperm, mmr = mmr, m = m))
}

#' Set random seed 
set.seed(124)

#' Generate correctly linked data (X)
n <- 150 # (!) reduced n to 150 for tractability (originally 1000)
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

### Unrelaxed JAGS 
library(rjags)

# (!) omit 'eta', 'alpha_prior', and 'm' as there is no grid approximation here.
dat <- list(
 x = X[,2], 
 z = X[,2], 
 y = yperm, 
 n = n
)

# Compile the exact model
jm <- jags.model(file = "ext1o-jags.txt", 
                 data = dat,
                 n.chains = 1, 
                 n.adapt = 1000, 
                 quiet = FALSE,
                 inits = list(beta0 = true_beta[1],
                              beta1 = true_beta[2]))

update(jm, 1000)

# (!) Track only 'sigma' since 'tau' belongs to the approximation model
samples <- coda.samples(jm, 
                        variable.names = c("beta0", "beta1", "gamma0", "gamma1", "sigma"), 
                        n.iter = 2000)

summary(samples)
plot(samples)
