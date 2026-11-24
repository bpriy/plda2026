library(dplyr)

# set working directory to survey_results1"
rds_files <- list.files(pattern = "^sim_result_task_\\d+\\.rds$")
if (length(rds_files) > 0) {
 results <- do.call(rbind, lapply(rds_files, readRDS))
}

if (!is.null(results) && nrow(results) > 0) {
 summary_list <- split(results, list(results$Method, results$Target_Rate))
 
 idx <- which(results$Method == "mix_unwt")
 
 # get the true dynamic gamma values used in each iteration
 true_g0 <- results$G0_Est[idx] - results$G0_Bias[idx]
 true_g1 <- results$G1_Est[idx] - results$G1_Bias[idx]
 
 # correct the estimates by flipping the sign
 results$G0_Est[idx] <- -results$G0_Est[idx]
 results$G1_Est[idx] <- -results$G1_Est[idx]
 
 # recalculate the true Bias
 results$G0_Bias[idx] <- results$G0_Est[idx] - true_g0
 results$G1_Bias[idx] <- results$G1_Est[idx] - true_g1
 
 # recalculate the true Confidence Interval Coverage
 results$G0_Cov[idx] <- as.numeric(abs(results$G0_Est[idx] - true_g0) <= 1.96 * results$G0_SE[idx])
 results$G1_Cov[idx] <- as.numeric(abs(results$G1_Est[idx] - true_g1) <= 1.96 * results$G1_SE[idx])
 
 summary_list <- split(results, list(results$Method, results$Target_Rate))
 
 final_tab <- do.call(rbind, lapply(summary_list, function(sub) {
  if(nrow(sub) == 0) return(NULL) 
  data.frame(
   Target_Rate = sub$Target_Rate[1], Method = sub$Method[1],
   B0_Bias = mean(sub$B0_Bias, na.rm=TRUE), 
   B0_SD = sd(sub$B0_Est, na.rm=TRUE), 
   B0_SE = mean(sub$B0_SE, na.rm=TRUE), 
   B0_Cov = mean(sub$B0_Cov, na.rm=TRUE),
   B1_Bias = mean(sub$B1_Bias, na.rm=TRUE), 
   B1_SD = sd(sub$B1_Est, na.rm=TRUE),
   B1_SE = mean(sub$B1_SE, na.rm=TRUE), 
   B1_Cov = mean(sub$B1_Cov, na.rm=TRUE),
   B2_Bias = mean(sub$B2_Bias, na.rm=TRUE), 
   B2_SD = sd(sub$B2_Est, na.rm=TRUE), 
   B2_SE = mean(sub$B2_SE, na.rm=TRUE), 
   B2_Cov = mean(sub$B2_Cov, na.rm=TRUE),
   B3_Bias = mean(sub$B3_Bias, na.rm=TRUE), 
   B3_SD = sd(sub$B3_Est, na.rm=TRUE), 
   B3_SE = mean(sub$B3_SE, na.rm=TRUE), 
   B3_Cov = mean(sub$B3_Cov, na.rm=TRUE),
   G0_Bias = mean(sub$G0_Bias, na.rm=TRUE), 
   G0_SD = sd(sub$G0_Est, na.rm=TRUE), 
   G0_SE = mean(sub$G0_SE, na.rm=TRUE), 
   G0_Cov = mean(sub$G0_Cov, na.rm=TRUE),
   G1_Bias = mean(sub$G1_Bias, na.rm=TRUE),
   G1_SD = sd(sub$G1_Est, na.rm=TRUE), 
   G1_SE = mean(sub$G1_SE, na.rm=TRUE), 
   G1_Cov = mean(sub$G1_Cov, na.rm=TRUE),
   Time = mean(sub$Time, na.rm=TRUE)
  )
 }))
 
 order_methods <- c("oracle", "naive_unwt", "naive_wt",
                    "naive_des", "mix_unwt", "mix_wt", 
                    "mix_des_sand", "mix_des_boot")
 final_tab <- final_tab[order(final_tab$Target_Rate, 
                              match(final_tab$Method, order_methods)), ]
 method_names <- c("oracle" = "Oracle", 
                   "naive_unwt" = "Naive Survey-Ignorant", 
                   "naive_wt" = "Naive Design-Ignorant", 
                   "naive_des" = "Naive Design-Adjusted", 
                   "mix_unwt" = "Mixture Survey-Ignorant", 
                   "mix_wt" = "Mixture Design-Ignorant", 
                   "mix_des_sand" = "Mixture Design-Adjusted (Sandwich)", 
                   "mix_des_boot" = "Mixture Design-Adjusted (Bootstrap)")
 final_tab$Method <- method_names[final_tab$Method]
 
 oracle_row <- final_tab[final_tab$Method == "Oracle", ][1, ]
 non_oracle_tab <- final_tab[final_tab$Method != "Oracle", ]
 fmt <- function(x) sprintf("%.3f", x)
 
 # beta0 and beta1
 cat("\\begin{table}[htbp]\n\\centering\n")
 cat("\\caption{Outcome Model Parameters ($\\beta_0$ and $\\beta_1$)}\n\\label{tab:sim_beta1}\n")
 cat("\\renewcommand{\\arraystretch}{1.3}\n")
 cat("\\resizebox{\\textwidth}{!}{\n")
 cat("\\begin{tabular}{|l|*{8}{c|}}\n\\hline\n")
 cat("\\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{4}{c|}{$\\widehat{\\beta}_0$} & \\multicolumn{4}{c|}{$\\widehat{\\beta}_1$} \\\\ \\cline{2-9}\n")
 cat(" & Bias & SD & SE & Cov & Bias & SD & SE & Cov \\\\ \\hline\n")
 
 cat("\\multicolumn{9}{|l|}{\\textbf{Baseline (No Linkage Error)}} \\\\ \\hline\n")
 row_str <- paste(fmt(as.numeric(oracle_row[c("B0_Bias", "B0_SD", "B0_SE", "B0_Cov", "B1_Bias", "B1_SD", "B1_SE", "B1_Cov")])), collapse = " & ")
 cat(sprintf("%s & %s \\\\ \\hline\n", oracle_row$Method, row_str))
 
 for (rate in unique(non_oracle_tab$Target_Rate)) {
  cat(sprintf("\\multicolumn{9}{|l|}{\\textbf{Mismatch Rate $\\alpha = %.2f$}} \\\\ \\hline\n", rate))
  sub_dat <- non_oracle_tab[non_oracle_tab$Target_Rate == rate, ]
  for (i in 1:nrow(sub_dat)) {
   row_str <- paste(fmt(as.numeric(sub_dat[i, c("B0_Bias", "B0_SD", "B0_SE", "B0_Cov", "B1_Bias", "B1_SD", "B1_SE", "B1_Cov")])), collapse = " & ")
   cat(sprintf("%s & %s \\\\ \\hline\n", sub_dat$Method[i], row_str))
  }
 }
 cat("\\end{tabular}\n}\n\\end{table}\n")
 
 # beta2 and beta3
 cat("\\begin{table}[htbp]\n\\centering\n")
 cat("\\caption{Outcome Model Parameters ($\\beta_2$ and $\\beta_3$)}\n\\label{tab:sim_beta2}\n")
 cat("\\renewcommand{\\arraystretch}{1.3}\n")
 cat("\\resizebox{\\textwidth}{!}{\n")
 cat("\\begin{tabular}{|l|*{8}{c|}}\n\\hline\n")
 cat("\\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{4}{c|}{$\\widehat{\\beta}_2$} & \\multicolumn{4}{c|}{$\\widehat{\\beta}_3$} \\\\ \\cline{2-9}\n")
 cat(" & Bias & SD & SE & Cov & Bias & SD & SE & Cov \\\\ \\hline\n")
 
 cat("\\multicolumn{9}{|l|}{\\textbf{Baseline (No Linkage Error)}} \\\\ \\hline\n")
 row_str <- paste(fmt(as.numeric(oracle_row[c("B2_Bias", "B2_SD", "B2_SE", "B2_Cov", "B3_Bias", "B3_SD", "B3_SE", "B3_Cov")])), collapse = " & ")
 cat(sprintf("%s & %s \\\\ \\hline\n", oracle_row$Method, row_str))
 
 for (rate in unique(non_oracle_tab$Target_Rate)) {
  cat(sprintf("\\multicolumn{9}{|l|}{\\textbf{Mismatch Rate $\\alpha = %.2f$}} \\\\ \\hline\n", rate))
  sub_dat <- non_oracle_tab[non_oracle_tab$Target_Rate == rate, ]
  for (i in 1:nrow(sub_dat)) {
   row_str <- paste(fmt(as.numeric(sub_dat[i, c("B2_Bias", "B2_SD", "B2_SE", "B2_Cov", "B3_Bias", "B3_SD", "B3_SE", "B3_Cov")])), collapse = " & ")
   cat(sprintf("%s & %s \\\\ \\hline\n", sub_dat$Method[i], row_str))
  }
 }
 cat("\\end{tabular}\n}\n\\end{table}\n")
 
 # gamma
 tab2_data <- final_tab[grepl("Mixture", final_tab$Method), ]
 cat("\\begin{table}[htbp]\n\\centering\n")
 cat("\\caption{Mismatch Indicator Model Parameters}\n\\label{tab:sim_gamma}\n")
 cat("\\renewcommand{\\arraystretch}{1.3}\n")
 cat("\\resizebox{\\textwidth}{!}{\n")
 cat("\\begin{tabular}{|l|*{8}{c|}}\n\\hline\n")
 cat("\\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{4}{c|}{$\\widehat{\\gamma}_0$ (Intercept)} & \\multicolumn{4}{c|}{$\\widehat{\\gamma}_1$ (Slope for Z)} \\\\ \\cline{2-9}\n")
 cat(" & Bias & SD & SE & Cov & Bias & SD & SE & Cov \\\\ \\hline\n")
 
 for (rate in unique(tab2_data$Target_Rate)) {
  cat(sprintf("\\multicolumn{9}{|l|}{\\textbf{Mismatch Rate $\\alpha = %.2f$}} \\\\ \\hline\n", rate))
  sub_dat <- tab2_data[tab2_data$Target_Rate == rate, ]
  for (i in 1:nrow(sub_dat)) {
   row_vals <- sub_dat[i, c("G0_Bias", "G0_SD", "G0_SE", "G0_Cov", "G1_Bias", "G1_SD", "G1_SE", "G1_Cov")]
   row_str <- paste(fmt(as.numeric(row_vals)), collapse = " & ")
   cat(sprintf("%s & %s \\\\ \\hline\n", sub_dat$Method[i], row_str))
  }
 }
 cat("\\end{tabular}\n}\n\\end{table}\n")
 
 # run time
 tab3 <- final_tab[, c("Target_Rate", "Method", "Time")]
 tab3_wide <- reshape(tab3, idvar = "Method", timevar = "Target_Rate", direction = "wide")
 tab3_wide <- tab3_wide[match(method_names, tab3_wide$Method), ] 
 
 cat("\\begin{table}[htbp]\n\\centering\n")
 cat("\\caption{Average Run Time (Seconds) Per Replication}\n\\label{tab:sim_time}\n")
 cat("\\renewcommand{\\arraystretch}{1.3}\n")
 cat("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}|l|c|c|c|}\n\\hline\n")
 cat(sprintf("\\textbf{Method} & \\textbf{$\\alpha = %.2f$} & \\textbf{$\\alpha = %.2f$} & \\textbf{$\\alpha = %.2f$} \\\\ \\hline\n", 
             unique(final_tab$Target_Rate)[1], unique(final_tab$Target_Rate)[2], unique(final_tab$Target_Rate)[3]))
 
 for (i in 1:nrow(tab3_wide)) {
  time_vals <- tab3_wide[i, 2:4]
  time_str <- paste(sprintf("%.2f", as.numeric(time_vals)), collapse = " & ")
  cat(sprintf("%s & %s \\\\ \\hline\n", tab3_wide$Method[i], time_str))
 }
 cat("\\end{tabular*}\n\\end{table}\n")
}