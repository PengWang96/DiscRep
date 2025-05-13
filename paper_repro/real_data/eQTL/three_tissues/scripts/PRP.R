rm(list = ls())
library(DiscRep)
library(ggplot2)
library(foreach)
library(doParallel)
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

t1 <- Sys.time()
set.seed(234)

data_MS <- read.table("../data/Muscle_Skeletal.torus.gene.eff", header = F); data_MS <- na.omit(data_MS)
data_AA <- read.table("../data/Artery_Aorta.torus.gene.eff", header = F); data_AA <- na.omit(data_AA)
data_L <- read.table("../data/Liver.torus.gene.eff", header = F); data_L <- na.omit(data_L)
head(data_MS)
head(data_AA)
head(data_L)

# 删除第三列或第四列中存在 0 的行
data_MS_filtered <- data_MS[data_MS$V3 != 0 & data_MS$V4 != 0, ]
data_AA_filtered <- data_AA[data_AA$V3 != 0 & data_AA$V4 != 0, ]
data_L_filtered <- data_L[data_L$V3 != 0 & data_L$V4 != 0, ]
# 通过 gene (V1) 和 snp (V2) 找到交集
common_genes <- Reduce(intersect, list(data_MS_filtered$V1, data_AA_filtered$V1, data_L_filtered$V1))
# 保留交集行
data_MS_final <- data_MS_filtered[data_MS_filtered$V1 %in% common_genes, ]
data_AA_final <- data_AA_filtered[data_AA_filtered$V1 %in% common_genes, ]
data_L_final <- data_L_filtered[data_L_filtered$V1 %in% common_genes, ]
data_MS_final <- data_MS_final[, -2]
data_AA_final <- data_AA_final[, -2]
data_L_final <- data_L_final[, -2]
colnames(data_MS_final) <- c("gene", "bhat_Muscle_Skeletal", "se_Muscle_Skeletal")
colnames(data_AA_final) <- c("gene", "bhat_Artery_Aorta", "se_Artery_Aorta")
colnames(data_L_final) <- c("gene", "bhat_Liver", "se_Liver")

# 按 gene 列合并数据框
data <- Reduce(function(x, y) merge(x, y, by = "gene"),
               list(data_AA_final, data_L_final, data_MS_final))
head(data)
write.table(data, file = "../data/simple_data.3tissue.summary",
            sep = "\t", row.names = FALSE, quote = FALSE)

subsample <- 1:nrow(data)
gene <- data[subsample, 1]
hat_beta <- as.matrix(data[subsample, c(2, 4, 6)])
hat_sigma_sq <- as.matrix(data[subsample, c(3, 5, 7)]^2)
rm(data_MS, data_AA, data_L,
   data_MS_filtered, data_AA_filtered, data_L_filtered,
   data_MS_final, data_AA_final, data_L_final, data)
gc()

pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate
m <- ncol(hat_beta)
num <- nrow(hat_beta)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

parts <- 1
rows_per_part <- floor(num / parts)
remainder <- num %% parts
results_list <- vector("list", parts)

for(i in 1:parts) {
  start_row <- ((i - 1) * rows_per_part) + 1
  end_row <- i * rows_per_part
  if (i == parts) {
    end_row <- end_row + remainder
  }

  part_hat_beta <- hat_beta[start_row:end_row, ]
  part_hat_sigma_sq <- hat_sigma_sq[start_row:end_row, ]

  sim_results <- foreach(x = 1:nrow(part_hat_beta), .combine = c, .packages = "DiscRep") %dopar% {
    metropolis_hastings(N, r, m, part_hat_beta[x,], part_hat_sigma_sq[x,], k_vec = k_vec)$p_value #
    # fixed_effect(N, m, part_hat_beta[x,], part_hat_sigma_sq[x,])$p_value # MC in FE
    # frequency_pvalue(part_hat_beta[x,], part_hat_sigma_sq[x,]) # Q test in FE
    # calc_PRP_FE_NumInt(part_hat_beta[x,], part_hat_sigma_sq[x,])$Integral_Value # numerically calc PRP in FE
  }
  results_list[[i]] <- sim_results
}

stopCluster(cl)
p_values <- unlist(results_list)


load("../output/results_PRP.rda")
p_values <- result$p_values
plot <-
ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.title.x = element_text(size = 22), # Increase x-axis title font size
    axis.title.y = element_text(size = 22)  # Increase y-axis title font size
  ) +
  xlab("Posterior-PRPs under Distinguishability Criterion") +
  ylab("Count")
plot
ggsave("../plot/eqtl_PRP_hist.pdf", plot, height = 8, width = 10)
dev.off()

result <- data.frame(gene = gene, p_values = p_values)
sum(result$p_values <= 0.05)
# save(result, file = "../output/results_PRP_Pmis0.05(MH-MCMC).rda")
# save(result, file = "../output/results_PRP_Pmis0(MH-MCMC).rda")
# save(result, file = "../output/results_PRP_Pmis0(MC).rda")
# save(result, file = "../output/results_Cochran_Q.rda")
# write.table(result, file = "../output/results_Cochran_Q.txt",
#             sep = "\t", row.names = FALSE, quote = FALSE)
save(result, file = "../output/results_PRP.rda") # without fixed heterogeneity_level
write.table(result, file = "../output/results_PRP.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
# save(result, file = "../output/results_PRP_Pmis0(NumInt).rda")
t2 <- Sys.time()
(t2 - t1)











# pvec = c(10^seq(-10, log10(0.05), 0.01),0.05)
# k_vec <- sapply(pvec, inverse_P_mis)
#
# pvec <- seq(0, 0.05, length = 1000)
# k_vec <- sapply(pvec, inverse_P_mis)
#
# rm(list = ls())
# library(DiscRep)
# library(ggplot2)
# library(foreach)
# library(doParallel)
# library(rstudioapi)
# library(rprojroot)
# if (interactive()) {
#   current_file <- rstudioapi::getActiveDocumentContext()$path
# } else {
#   args <- commandArgs(trailingOnly = FALSE)
#   current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
# }
# setwd(dirname(current_file))
# print(getwd())
#
# t1 <- Sys.time()
# set.seed(234)
# data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
# head(data)
#
#
# subsample <- 1:nrow(data)
# # subsample <- which(data$Gene == "ENSG00000176945.16")
# gene <- data[subsample, 1]
# hat_beta <- as.matrix(data[subsample, c(2, 4, 6)])
# hat_sigma_sq <- as.matrix(data[subsample, c(3, 5, 7)]^2)
#
# P_mis_vec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
# k_vec <- sapply(P_mis_vec, inverse_P_mis)
# gc()
#
# N <- 10000
# r <- 0.05 # burn-in rate
# m <- ncol(hat_beta)
# num <- nrow(hat_beta)
#
# numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
#
# parts <- 1
# rows_per_part <- floor(num / parts)
# remainder <- num %% parts
# results_list <- vector("list", length(P_mis_vec))  # 保存每个P_mis_vec对应的结果
#
#
# for(p_idx in 1:length(P_mis_vec)) {
#   heterogeneity_level <- k_vec[p_idx]
#   p_values_per_pmis <- vector("list", parts)  # 每个parts的结果将存储在这个list中
#
#   for(i in 1:parts) {
#     start_row <- ((i - 1) * rows_per_part) + 1
#     end_row <- i * rows_per_part
#     if (i == parts) {
#       end_row <- end_row + remainder
#     }
#     part_hat_beta <- hat_beta[start_row:end_row, ]
#     part_hat_sigma_sq <- hat_sigma_sq[start_row:end_row, ]
#
#     p_values_per_pmis[[i]] <- foreach(x = 1:nrow(part_hat_beta), .combine = c, .packages = "DiscRep") %dopar% {
#       metropolis_hastings(N, r, m, part_hat_beta[x, ], part_hat_sigma_sq[x, ], heterogeneity_level = heterogeneity_level)$p_value
#     }
#   }
#   results_list[[p_idx]] <- do.call(c, p_values_per_pmis)
# }
# stopCluster(cl)
#
# p_values_matrix <- do.call(cbind, results_list)
# colnames(p_values_matrix) <- paste0("P_mis_", as.character(P_mis_vec))
# p_values_df <- as.data.frame(p_values_matrix)
# head(p_values_df)
# #   xlab("Frequentist p-values") +
#
#
# result <- as.data.frame(cbind(gene, p_values_df))
# save(result, file = "../output/results_PRP_matrix.rda") # without fixed heterogeneity_level
# write.table(result, file = "../output/results_PRP_matrix.txt",
#             sep = "\t", row.names = FALSE, quote = FALSE)
# t2 <- Sys.time()
# (t2 - t1)
# ggplot(result, aes(x = P_mis_0.01)) +
#   geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
#   scale_x_continuous(limits = c(0, 1)) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 12, hjust = 0.5)) +
#   xlab("Posterior p-values under distinguishability criterion") +
#   ylab("Count")
