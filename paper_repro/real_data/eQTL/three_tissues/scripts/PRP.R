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
    metropolis_hastings(N, r, m, part_hat_beta[x,], part_hat_sigma_sq[x,], k_vec = k_vec)$p_value
  }
  results_list[[i]] <- sim_results
}

stopCluster(cl)
p_values <- unlist(results_list)

result <- data.frame(gene = gene, p_values = p_values)
sum(result$p_values <= 0.05)
saveRDS(result, file = "../output/results_PRP.rds") # without fixed k
write.table(result, file = "../output/results_PRP.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
t2 <- Sys.time()
(t2 - t1)