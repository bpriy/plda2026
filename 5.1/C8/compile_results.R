
library(purrr)
library(dplyr)

# Summary Functions
bias_fun <- function(val, theta) { colMeans(val[,1:length(theta)], na.rm=TRUE) - theta }
sd_fun <- function(val) { apply(val, 2, sd, na.rm=TRUE) }
cov_fun <- function(thetahat, se, theta){
 cp <- numeric(length(theta))
 for (i in 1:length(theta)){
  lower <- thetahat[,i] - qt(0.975, df = 995)*se[,i]
  upper <- thetahat[,i] + qt(0.975, df = 995)*se[,i]
  cov <- theta[i] >= lower & theta[i] <= upper
  cp[i] <- mean(cov, na.rm = TRUE)
 }
 return(cp)
}

SETTINGS <- c("Step1_Null", "Step2_PosCorr", 
              "Step3_K5", "Step3_K10", "Step3_K20", 
              "Step4_E0.1", "Step4_E0.5", "Step4_E1.0", 
              "Step5A_Probit", "Step5B_Cloglog", "Step5C_Interact")

Truebeta <- c(1, 2, -1.5, 1)
Truesigma <- 0.25

# Compile all results
for (setting in SETTINGS) {
 cat("RESULTS FOR SETTING:", setting, "\n")
 
 # Adjust target parameters based on the setting
 Truegamma <- c(1, -2.5, 1)
 if (setting == "Step1_Null") Truegamma <- c(1, 0, 1)
 theta <- c(Truebeta, Truesigma, Truegamma)
 
 out_dir <- paste0("results_", setting)
 if(!dir.exists(out_dir)) { cat("Directory missing.\n"); next }
 
 methods <- c("naive", "oracle", "aint", "adelta", "az", "ext1", "ext1z", "aoracle")
 
 bias_res <- matrix(NA, nrow=length(methods), ncol=length(theta))
 sd_res <- matrix(NA, nrow=length(methods), ncol=length(theta))
 cr_res <- matrix(NA, nrow=length(methods), ncol=length(theta))
 time_res <- matrix(NA, nrow=length(methods), ncol=1)
 
 rownames(bias_res) <- rownames(sd_res) <- rownames(cr_res) <- rownames(time_res) <- methods
 colnames(bias_res) <- colnames(sd_res) <- colnames(cr_res) <- c("b0", "b1", "b2", "b3", "sigma", "g0", "g1", "g2")
 colnames(time_res) <- c("Avg_Time_Secs")
 
 for (m in methods) {
  # Read coefficients
  coef_files <- list.files(out_dir, pattern = paste0("^", m, "_coef"), full.names = TRUE)
  if(length(coef_files) == 0) next
  coef_val <- do.call(rbind, map(coef_files, readRDS))
  
  # Read standard errors
  se_files <- list.files(out_dir, pattern = paste0("^", m, "_se"), full.names = TRUE)
  se_val <- do.call(rbind, map(se_files, readRDS))
  
  # Read run times
  time_files <- list.files(out_dir, pattern = paste0("^", m, "_time"), full.names = TRUE)
  if(length(time_files) > 0) {
   time_val <- do.call(rbind, map(time_files, readRDS))
   time_res[m, 1] <- mean(time_val, na.rm=TRUE)
  }
  
  # Delta method correction for Plain adjustments
  if(m %in% c("aint", "adelta", "az", "aoracle")){
   se_val[,5] <- se_val[,5]/(2*coef_val[,5]) 
  }
  
  bias_res[m, ] <- bias_fun(coef_val, theta)
  sd_res[m, ]   <- sd_fun(coef_val)
  cr_res[m, ]   <- cov_fun(coef_val, se_val, theta)
 }
 
 cat("\n BIAS \n"); print(round(bias_res, 4))
 cat("\n EMPIRICAL SD \n"); print(round(sd_res, 4))
 cat("\n 95% COVERAGE RATE (CR) \n"); print(round(cr_res, 3))
 cat("\n AVG RUN TIME (s) \n"); print(round(time_res, 4))
}