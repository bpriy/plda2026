# Load required libraries
library(dplyr)
library(ggplot2)
library(scales)
library(purrr)

results_dir <- "results_scalability/N_scaling"

# Read and combine all batch .rds files
rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(rds_files, readRDS)

# Summarize the empirical runtimes across the 1,000 replications
# We calculate the mean runtime and the 95% Confidence Interval bounds
scaling_summary <- all_results %>%
 group_by(N) %>%
 summarize(
  Mean_Runtime = mean(Runtime_Secs),
  SD_Runtime = sd(Runtime_Secs),
  Count = n(),
  SE_Runtime = SD_Runtime / sqrt(Count),
  CI_Lower = Mean_Runtime - 1.96 * SE_Runtime,
  CI_Upper = Mean_Runtime + 1.96 * SE_Runtime,
  .groups = "drop"
 )

# Generate the Theoretical O(n) Reference Line
# Anchor the theoretical line to the empirical mean of the smallest sample size (N = 1000)
anchor_N <- 1000
anchor_time <- scaling_summary$Mean_Runtime[scaling_summary$N == anchor_N]

# Create a data frame for the theoretical line (Slope = 1 on Log-Log scale)
theoretical_line <- data.frame(
 N = scaling_summary$N,
 Theoretical_Time = anchor_time * (scaling_summary$N / anchor_N)
)

# Build the Plot
pub_plot <- ggplot() +
 geom_line(data = theoretical_line, 
           aes(x = N, y = Theoretical_Time, color = "Theoretical O(n)", linetype = "Theoretical O(n)"), 
           linewidth = 1) +
 geom_line(data = scaling_summary, 
           aes(x = N, y = Mean_Runtime, color = "Empirical Runtime", linetype = "Empirical Runtime"), 
           linewidth = 1.2) +
 geom_errorbar(data = scaling_summary, 
               aes(x = N, 
                   ymin = Mean_Runtime - SD_Runtime, 
                   ymax = Mean_Runtime + SD_Runtime, 
                   color = "Empirical Runtime"), 
               width = 0.15, linewidth = 1.2) +
 geom_point(data = scaling_summary, 
            aes(x = N, y = Mean_Runtime, color = "Empirical Runtime"), 
            size = 2.5) +
 scale_x_log10(
  breaks = c(10^3, 10^4, 10^5, 10^6),
  labels = trans_format("log10", math_format(10^.x)),
  name = "Sample Size (n)"
 ) +
 scale_y_log10(
  breaks = 10^(-1:3),
  labels = trans_format("log10", math_format(10^.x)),
  limits = c(0.2, 1000), 
  name = "Runtime (Seconds)"
 ) +
 annotation_logticks(sides = "l") + 
 
 scale_color_manual(name = "", values = c("Empirical Runtime" = "blue", "Theoretical O(n)" = "gray40")) +
 scale_linetype_manual(name = "", values = c("Empirical Runtime" = "dashed", "Theoretical O(n)" = "solid")) +
 
 guides(color = guide_legend(override.aes = list(
  shape = c(16, NA),
  size = c(5, NA),
  linetype = c("dashed", "solid")
 )), linetype = "none") +
 theme_classic(base_size = 20) +
 theme(
  legend.position = c(0.25, 0.85),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.key = element_blank(),
  legend.key.width = unit(2.5, "cm"),
  legend.text = element_text(size = 16),
  axis.title = element_text(face = "bold", size = 20),
  axis.text = element_text(color = "black", size = 18),
  axis.ticks.length = unit(0.2, "cm")
 )

ggsave("KDE_runtimes_plot.pdf", plot = pub_plot, width = 7, height = 5, dpi = 300)
ggsave("KDE_runtimes_plot.png", plot = pub_plot, width = 7, height = 5, dpi = 300)

print(pub_plot)

print(scaling_summary)

# Parameter Evaluation (Internal Validation)

# Define True Parameters from the Simulation DGP
true_b0 <- 1
true_b1 <- -1
true_g0 <- 2.5
true_g1 <- 4.5
true_sig <- 0.25

# Compute Bias, SD, Mean SE, and CR for each parameter grouped by N
parameter_evaluation <- all_results %>%
 group_by(N) %>%
 summarize(
  B0_Bias = mean(Beta_0) - true_b0,
  B0_SD = sd(Beta_0),
  B0_Mean_SE = mean(SE_Beta_0),
  B0_CR = mean((Beta_0 - 1.96 * SE_Beta_0 <= true_b0) & 
                (Beta_0 + 1.96 * SE_Beta_0 >= true_b0)),
  
  B1_Bias = mean(Beta_1) - true_b1,
  B1_SD = sd(Beta_1),
  B1_Mean_SE = mean(SE_Beta_1),
  B1_CR = mean((Beta_1 - 1.96 * SE_Beta_1 <= true_b1) & 
                (Beta_1 + 1.96 * SE_Beta_1 >= true_b1)),
  
  G0_Bias = mean(Gamma_0) - true_g0,
  G0_SD = sd(Gamma_0),
  G0_Mean_SE = mean(SE_Gamma_0),
  G0_CR = mean((Gamma_0 - 1.96 * SE_Gamma_0 <= true_g0) & 
                (Gamma_0 + 1.96 * SE_Gamma_0 >= true_g0)),
  
  G1_Bias = mean(Gamma_1) - true_g1,
  G1_SD = sd(Gamma_1),
  G1_Mean_SE = mean(SE_Gamma_1),
  G1_CR = mean((Gamma_1 - 1.96 * SE_Gamma_1 <= true_g1) & 
                (Gamma_1 + 1.96 * SE_Gamma_1 >= true_g1)),
  
  Sig_Bias = mean(Sigma) - true_sig,
  Sig_SD = sd(Sigma),
  Sig_Mean_SE = mean(SE_Sigma),
  Sig_CR = mean((Sigma - 1.96 * SE_Sigma <= true_sig) & 
                 (Sigma + 1.96 * SE_Sigma >= true_sig)),
  
  .groups = "drop"
 )

cat("\n Parameter Evaluation (Bias, SD, SE, Coverage) \n")
print(as.data.frame(t(parameter_evaluation)))