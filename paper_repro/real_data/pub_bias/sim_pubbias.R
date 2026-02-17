rm(list = ls())
t1 <- Sys.time()
library(metafor)
library(DiscRep)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())


# Save funnel plots for all three datasets
datasets <- list(
  # Stead et al. (2012)
  list(name = "dat.slf", file = "./plot/funnel_slf.pdf"),
  # Hrobjartsson and Gotzsche (2010)
  list(name = "dat.ha", file = "./plot/funnel_ha.pdf"),
  # Liu CJ, Latham NK (2009)
  list(name = "dat.lcj", file = "./plot/funnel_lcj.pdf")
)

pdf_w <- 7
pdf_h <- 5

for (ds in datasets) {
  data(list = ds$name)
  data <- get(ds$name)
  res <- rma(yi = y, sei = sqrt(s2), data = data)
  cairo_pdf(ds$file, width = pdf_w, height = pdf_h)
  par(mar = c(4.5, 5.5, 2, 2) + 0.1)
  par(cex.axis = 2, cex.lab = 2.2)
  funnel(res)
  dev.off()
  message("Saved: ", ds$file)
}


# Run analysis for all three datasets
pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
k_vec <- sapply(pvec, inverse_P_mis)

for (ds in datasets) {
  cat("\n========== ", ds$name, " ==========\n")
  data <- get(ds$name)
  res <- rma(yi = y, sei = sqrt(s2), data = data)
  print(res)

  # Egger's test
  egger_test <- regtest(res, model = "lm", predictor = "sei")
  cat("Egger test:\n")
  print(egger_test)

  # Begg's test
  begg_test <- ranktest(res)
  cat("Begg test:\n")
  print(begg_test)

  # Metropolis-Hastings
  m <- nrow(data)
  hat_beta <- data$y
  hat_sigma_sq <- data$s2
  sim_results_q <- metropolis_hastings(10000, 0.05, m,
    hat_beta, hat_sigma_sq,
    test = "Q", k_vec = k_vec
  )
  cat("MH Q test p-value:", sim_results_q$p_value, "\n")
  sim_results_egger <- metropolis_hastings(10000, 0.05, m,
    hat_beta, hat_sigma_sq,
    test = "Egger", k_vec = k_vec
  )
  cat("MH Egger test p-value:", sim_results_egger$p_value, "\n")

  # Q test
  q_pval <- frequency_pvalue(hat_beta, hat_sigma_sq)
  cat("Q test p-value:", q_pval, "\n")
}
t2 <- Sys.time()
(t2 - t1)
