mismatches_fun <- function(gammaval){
 yperm <- y
 
 #' generate mismatches (w/ z = X)
 m <- 1 - rbinom(n, prob = plogis(X %*% gammaval), size = 1) 
 #m <- rbinom(n, size = 1, prob = 0.3)
 mmr <- sum(m)/n
 
 yperm <- y
 shuffled_ix <- sample(which(m == 1)) # only shuffle the mismatched data
 yperm[shuffled_ix] <- yperm[c(shuffled_ix[2:length(shuffled_ix)], 
                               shuffled_ix[1])]
 
 return(list(yperm = yperm, mmr = mmr, m = m))
}

set.seed(124)

#' Generate correctly linked data (X)
n <- 1000
p <- 1 # number of covariates
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X <- cbind(1, X) # intercept + covariates

beta <- c(1, -1)
gamma <- c(2.5, 4.5) 
sigma <- 0.25

y <- rnorm(n, X %*% beta, sigma)
mmdata <- mismatches_fun(gamma)

m <- mmdata$m
g <- function(yval){dnorm(yval, mean = mean(y[m == 0]), sd = sd(y[m == 0]))} # Gaussian density
g1 <- function(yval){dnorm(yval, mean = mean(y[m == 1]), sd = sd(y[m == 1]))} # Gaussian density
g2 <- function(yval){dnorm(yval, mean = mean(y), sd = sd(y))} # Gaussian density

par(mfrow = c(1,3))
plot(X[,2], mmdata$yperm, col = "red", pch = 17, xlab = "Predictor x", 
     ylab = "Response y", cex = 1.5, cex.lab=1.2)
points(X[,2], y, pch = 20, cex = 1.5)

plot(y, mmdata$yperm, col = "black", xlab = "True y", ylab = "Permuted y", 
     pch = 20, cex = 1.5, cex.lab=1.2)
abline(a = 0, b= 1, col = "red", lwd = 4, lty = "dashed", cex = 1.4)

curve(g1, -2, 5, xlab = "y", ylab = "fy", cex.lab=1.2, lwd = 2, lty = "dashed")
curve(g, -2, 5, add = TRUE, lwd = 2)
curve(g2, -2, 5, add = TRUE, type = "b")
# 10 x 3