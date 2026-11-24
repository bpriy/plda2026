library(dplyr)
library(kableExtra)

# Table 1: 6 Methodologies
m_files <- list.files("results_scalability/methods", pattern = "\\.rds$", full.names = TRUE)
df_methods <- do.call(rbind, lapply(m_files, readRDS))

tbl1_data <- df_methods %>%
 group_by(Method) %>%
 summarize(
  Beta_0 = sprintf("%.4f \n (%.4f)", mean(Beta_0), mean(SE_Beta_0)),
  Beta_1 = sprintf("%.4f \n (%.4f)", mean(Beta_1), mean(SE_Beta_1)),
  Sigma = sprintf("%.4f \n (%.4f)", mean(Sigma), mean(SE_Sigma)),
  Gamma_0 = sprintf("%.4f \n (%.4f)", mean(Gamma_0), mean(SE_Gamma_0)),
  Gamma_1 = sprintf("%.4f \n (%.4f)", mean(Gamma_1), mean(SE_Gamma_1)),
  `Runtime (s)` = sprintf("%.4f", mean(Runtime)),
  `EM Iters` = round(mean(Iters)),
  `Monotonic?` = ifelse(all(Monotonic), "Yes", "No")
 )

true_row <- data.frame(
 Method = "True Data Generating Parameters",
 Beta_0 = "1", Beta_1 = "-1", Sigma = "0.25", Gamma_0 = "2.5", Gamma_1 = "4.5",
 `Runtime (s)` = "NA", `EM Iters` = "NA", `Monotonic?` = "NA",
 check.names = FALSE
)
tbl1_data <- rbind(true_row, tbl1_data)

order_levels <- c("True Data Generating Parameters", "Original Dense O(n^2)", 
                  "Sparse (Thresholding)", "Block-Diagonal", "Exact O(n) Analytical", 
                  "Exact + KDE", "Exact + KDE + Original")
tbl1_data <- tbl1_data[match(order_levels, tbl1_data$Method), ]

colnames(tbl1_data)[1:6] <- c("Methodology", "$\\hat{\\beta}_0$ (SE)", 
                              "$\\hat{\\beta}_1$ (SE)", 
                              "$\\hat{\\sigma}$ (SE)", 
                              "$\\hat{\\gamma}_0$ (SE)", 
                              "$\\hat{\\gamma}_1$ (SE)")
 
tbl1_data %>%
 kbl(format = "html", escape = FALSE, align = c("l", rep("c", 8))) %>%
 kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
 row_spec(1, bold = TRUE, background = "#f2f2f2") %>%
 add_header_above(c(" " = 1, "Outcome Model" = 4, "Linkage Model" = 2, "Performance Diagnostics" = 3)) %>%
 print()

# Table 2: Blocks Correctly Specified
c_files <- list.files("results_scalability/correct", pattern = "\\.rds$", full.names = TRUE)
df_correct <- do.call(rbind, lapply(c_files, readRDS))

tbl2_data <- df_correct %>%
 group_by(Blocks) %>%
 summarize(
  Pairs_Evaluated = sprintf("%.0e", mean(Pairs_Evaluated)),
  Beta_0 = sprintf("%.4f \n (%.4f)", mean(Beta_0), mean(SE_Beta_0)),
  Beta_1 = sprintf("%.4f \n (%.4f)", mean(Beta_1), mean(SE_Beta_1)),
  Sigma = sprintf("%.4f \n (%.4f)", mean(Sigma), mean(SE_Sigma)),
  Gamma_0 = sprintf("%.4f \n (%.4f)", mean(Gamma_0), mean(SE_Gamma_0)),
  Gamma_1 = sprintf("%.4f \n (%.4f)", mean(Gamma_1), mean(SE_Gamma_1)),
  Runtime_Secs = sprintf("%.4f", mean(Runtime_Secs)),
  Iterations = round(mean(Iterations))
 )

colnames(tbl2_data)[3:7] <- c("$\\hat{\\beta}_0$ (SE)", "$\\hat{\\beta}_1$ (SE)", "$\\hat{\\sigma}$ (SE)", "$\\hat{\\gamma}_0$ (SE)", "$\\hat{\\gamma}_1$ (SE)")

tbl2_data %>%
 kbl(format = "html", escape = FALSE, align = c("l", rep("c", 8))) %>%
 kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
 print()

# Table 3: Blocks Misspecified
ms_files <- list.files("results_scalability/misspec", pattern = "\\.rds$", full.names = TRUE)
df_misspec <- do.call(rbind, lapply(ms_files, readRDS))

tbl3_data <- df_misspec %>%
 group_by(Assumed_Blocks) %>%
 summarize(
  Pairs_Evaluated = sprintf("%.0e", mean(Pairs_Evaluated)),
  Beta_0 = sprintf("%.4f \n (%.4f)", mean(Beta_0), mean(SE_Beta_0)),
  Beta_1 = sprintf("%.4f \n (%.4f)", mean(Beta_1), mean(SE_Beta_1)),
  Sigma = sprintf("%.4f \n (%.4f)", mean(Sigma), mean(SE_Sigma)),
  Gamma_0 = sprintf("%.4f \n (%.4f)", mean(Gamma_0), mean(SE_Gamma_0)),
  Gamma_1 = sprintf("%.4f \n (%.4f)", mean(Gamma_1), mean(SE_Gamma_1)),
  Time_Secs = sprintf("%.4f", mean(Time_Secs)),
  Iter = round(mean(Iter))
 )

colnames(tbl3_data)[3:7] <- c("$\\hat{\\beta}_0$ (SE)", "$\\hat{\\beta}_1$ (SE)", "$\\hat{\\sigma}$ (SE)", "$\\hat{\\gamma}_0$ (SE)", "$\\hat{\\gamma}_1$ (SE)")

tbl3_data %>%
 kbl(format = "html", escape = FALSE, align = c("l", rep("c", 8))) %>%
 kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
 print()