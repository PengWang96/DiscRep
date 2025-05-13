################################################################################
#######             Prepare data for eqtlbma            ########################
################################################################################
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

data_AA_final$snp <- paste0("snp", 1:nrow(data_AA_final))
data_AA_final$sigmahat <- 1
data_MS_final$snp <- paste0("snp", 1:nrow(data_MS_final))
data_MS_final$sigmahat <- 1
data_L_final$snp <- paste0("snp", 1:nrow(data_L_final))
data_L_final$sigmahat <- 1
identical(data_AA_final$V1, data_MS_final$V1)
identical(data_AA_final$V1, data_L_final$V1)

colnames(data_MS_final) <- c("gene", "n",	"betahat.geno",	"sebetahat.geno", "snp", "sigmahat")
colnames(data_AA_final) <- c("gene", "n",	"betahat.geno",	"sebetahat.geno", "snp", "sigmahat")
colnames(data_L_final) <- c("gene", "n",	"betahat.geno",	"sebetahat.geno", "snp", "sigmahat")


write.table(data_MS_final,
            file = gzfile("../eqtlbma/sstats_muscle_skeletal.txt.gz"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(data_AA_final,
            file = gzfile("../eqtlbma/sstats_artery_aorta.txt.gz"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(data_L_final,
            file = gzfile("../eqtlbma/sstats_liver.txt.gz"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
t2 <- Sys.time()
(t2 - t1)






################################################################################
#######       Compare p-values with PRP and eqtlbma     ########################
################################################################################
rm(list = ls())
library(ggplot2)
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
load("../output/results_PRP.rda")
PRP <- result
# load("../output/results_Cochran_Q.rda")
# Q_pval <- result
# load("../output/results_PRP_Pmis0(MH-MCMC).rda")
# PRP_k0 <- result
eqtlbma <- read.table("../eqtlbma/marginal_probability.txt.gz", header=T)

index <- match(PRP[, 1], eqtlbma[, 1])
eqtlbma <- eqtlbma[index, ]

if(identical(PRP[, 1], eqtlbma[, 1])) {
  print("The first column of a and b are identical.")
} else {
  print("The first column of a and b are not identical.")
}
# eqtlbma[ , 3:ncol(eqtlbma)] <- lapply(eqtlbma[ , 3:ncol(eqtlbma)], as.numeric)

eqtlbmaa <- 1 - eqtlbma$gene.post + eqtlbma$marg.config.1.2.3
# eqtlbmaa <- 1 - eqtlbmaa


PRPP <- PRP$p_values
2 * mean(PRPP)
2 * mean(eqtlbmaa)
sum(PRPP <= 0.05)/length(PRPP)
sum(eqtlbmaa <= 0.05)/length(eqtlbmaa)

cor.test(PRPP, eqtlbmaa)
plot(PRPP, eqtlbmaa)
plot(-log10(PRPP), -log10(eqtlbmaa))
cor.test(-log10(PRPP), -log10(eqtlbmaa))



library(RColorBrewer)
# my_colramp <- colorRampPalette(c("white", "#1f77b4"))
my_colramp <- colorRampPalette(brewer.pal(9, "YlGnBu"))
# my_colramp <- colorRampPalette(brewer.pal(9, "GnBu"))
# 设置全局字体大小
par(cex.lab = 2.2,     # x 和 y 轴标签字体大小
    cex.axis = 1.4,    # 坐标轴刻度字体大小
    cex.main = 1.8,    # 主标题字体大小
    mar = c(5, 5, 2, 2)) # 调整边距 (bottom, left, top, right)

smoothScatter(
  PRPP,
  eqtlbmaa,
  colramp = my_colramp,
  nrpoints = 0,
  transformation = sqrt,
  xlab = "Posterior-PRPs",
  ylab = "Tissue-consistency Probability (eQtlBma)"
)





data_comp <- data.frame(
  value = c(PRPP, eqtlbmaa),
  Method = factor(c(rep("posterior-PRP", length(PRPP)), rep("eqtlbma", length(eqtlbmaa))),
                    levels = c("posterior-PRP", "eqtlbma"))
)
# data_comp <- data.frame(
#   value = c(eqtlbma, PRP),
#   Method = factor(c(rep("eqtlbma", length(eqtlbma)), rep("posterior-PRP", length(PRP))),
#                   levels = c("eqtlbma", "posterior-PRP"))
# )

ggplot(data_comp, aes(x = value, fill = Method)) +
  geom_histogram(position = "dodge", bins = 20, alpha = 0.7, color = "black") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) + # 可自定义颜色
  labs(title = "Histogram Comparison of PRPP and eqtlbmaa",
       x = "Value",
       y = "Frequency",
       fill = "Method") +
  theme_bw()

ggplot(data.frame(eqtlbmaa), aes(x = eqtlbmaa)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.00, 1.00))

ggplot(PRP, aes(x = p_values)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.00, 1.00))




################################################################################
#######                     Forest plot                 ########################
################################################################################
library(forestplot)
library(ggplotify)
library(patchwork)
data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
data <- na.omit(data)
head(data)
# selected_genes <- data$gene[which(eqtlbmaa > 0.99 & PRPP < 0.001)] # [c(5, 8, 9, 11, 12)]
# selected_genes <- data$gene[which(eqtlbmaa > 0.99 & PRPP < 0.001)[c(5, 8, 9, 11, 12)]] #
selected_genes <- data$gene[which(eqtlbmaa < 0.001 & PRPP > 0.7)] #
selected_genes
length(selected_genes)
gene_symbols <- c("NBPF3", "GATD3", "PSPHP1", "MAP3K2-DT")
# Specify the number of plots you want (e.g., 3)
num_plots <- 4

# Generate the plots
plots <- list()
for (i in 1:num_plots) {
  selected_gene <- data[data$gene == selected_genes[i],]
  plot_matrix <- matrix(c(
    selected_gene$bhat_Artery_Aorta, selected_gene$bhat_Liver, selected_gene$bhat_Muscle_Skeletal,
    selected_gene$bhat_Artery_Aorta - 1.96 * selected_gene$se_Artery_Aorta,
    selected_gene$bhat_Liver - 1.96 * selected_gene$se_Liver,
    selected_gene$bhat_Muscle_Skeletal - 1.96 * selected_gene$se_Muscle_Skeletal,
    selected_gene$bhat_Artery_Aorta + 1.96 * selected_gene$se_Artery_Aorta,
    selected_gene$bhat_Liver + 1.96 * selected_gene$se_Liver,
    selected_gene$bhat_Muscle_Skeletal + 1.96 * selected_gene$se_Muscle_Skeletal
  ), nrow=3, byrow=TRUE)

  labels <- c('Artery Aorta', 'Liver', 'Muscle Skeletal')

  p <- forestplot(labeltext=labels,
                  mean=plot_matrix[1,],
                  lower=plot_matrix[2,],
                  upper=plot_matrix[3,],
                  xlab=paste("Effect Size for Gene", gene_symbols[i]),
                  lwd.ci = 3, # Make the lines bolder
                  txt_gp = fpTxtGp(
                    label = gpar(cex = 1.6),    # 调整y轴标签字体大小
                    ticks = gpar(cex = 1.2),    # 调整x轴刻度标签字体大小
                    xlab = gpar(cex = 1.6)      # 调整x轴标签字体大小
                  ),
                  boxsize = 0.12,               # 调整图形元素大小
                  lineheight = unit(2.2, "cm"), # 增加行高
                  graphwidth = unit(80, "mm"),  # 增加图形部分的宽度
                  # col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
                  new_page = F)
  plots[[i]] <- p
}

# Convert list of plots into grobs
plot_grobs <- lapply(plots, function(p) grid2grob(print(p)))

# Combine plots dynamically
p_combined <- wrap_elements(plot_grobs[[1]])
for (i in 2:length(plot_grobs)) {
  p_combined <- p_combined + wrap_elements(plot_grobs[[i]])
}
p_combined



#################### single gene ##########################
rm(list = ls())
library(ggplot2)
library(forestplot)
library(ggplotify)
library(patchwork)
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

selected_genes <- c(
  "ENSG00000112877.7",
  # "ENSG00000057468.6",
  "ENSG00000123200.16"
)
gene_symbols <- c(
  "CEP72",
  # "MSH4",
  "ZC3H13"
  )
gene_index <- 1

data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
data <- na.omit(data)
head(data)
selected_gene <- data[data$gene == selected_genes[gene_index],]
plot_matrix <- matrix(c(
  selected_gene$bhat_Artery_Aorta, selected_gene$bhat_Liver, selected_gene$bhat_Muscle_Skeletal,
  selected_gene$bhat_Artery_Aorta - 1.96 * selected_gene$se_Artery_Aorta,
  selected_gene$bhat_Liver - 1.96 * selected_gene$se_Liver,
  selected_gene$bhat_Muscle_Skeletal - 1.96 * selected_gene$se_Muscle_Skeletal,
  selected_gene$bhat_Artery_Aorta + 1.96 * selected_gene$se_Artery_Aorta,
  selected_gene$bhat_Liver + 1.96 * selected_gene$se_Liver,
  selected_gene$bhat_Muscle_Skeletal + 1.96 * selected_gene$se_Muscle_Skeletal
), nrow=3, byrow=TRUE)
labels <- c('Artery Aorta', 'Liver', 'Muscle Skeletal')

x_label <- as.expression(bquote("Effect Size for Gene" ~ italic(.(gene_symbols[gene_index]))))
clip <- list(c(-1.8, 1.8), c(-0.5, 0.5))
xticks <- list(c(-1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8),
               seq(-0.5, 0.5, by = 0.1))

p <- forestplot(labeltext = labels,
                mean = plot_matrix[1,],
                lower = plot_matrix[2,],
                upper = plot_matrix[3,],
                xlab = x_label,
                clip = clip[[gene_index]],        # 设置x轴范围
                xticks = xticks[[gene_index]],  # 自动生成刻度
                lwd.ci = 5,
                txt_gp = fpTxtGp(
                  label = gpar(cex = 1.9),
                  ticks = gpar(cex = 1.6),
                  xlab = gpar(cex = 1.9)
                ),
                boxsize = 0.1,
                lineheight = unit(6, "cm"),
                graphwidth = unit(200, "mm"),
                new_page = TRUE,
                col = fpColors(
                  box = "black",
                  line = "gray50",
                  zero = "gray50"
                ))
p





################################################################################
####### Compare varying Pmis ###################################################
################################################################################
rm(list = ls())
library(ggplot2)
setwd("D:/R/Replicability/Peng/real data/eQTL/three tissues/output")
load_data <- function(k) {
  load(paste0("reslutQ_k", k, ".rda"))
  data.frame(k = k, p_value = result_all$p_values)
}
k_values <- c(0, 0.2726814, 0.3349839, 0.4626489, 0.5420572, 100)
P_mis_values <- c(0, 0.05, 0.1, 0.2, 0.25, 0.5)
data <- do.call(rbind, lapply(k_values, load_data))

library(dplyr)
groups <- c("P_mis = 0, # of sig. p-value = 8439",
            "P_mis = 0.05, # of sig. p-value = 8286",
            "P_mis = 0.1, # of sig. p-value = 8053",
            "P_mis = 0.2, # of sig. p-value = 7412",
            "P_mis = 0.25, # of sig. p-value = 6913",
            "P_mis = 0.5, # of sig. p-value = 0")
data <- data %>%
  mutate(group = factor(k, levels = k_values, labels = groups))



# # Define a function to calculate the count of p-values less than 0.05
# count_p_less_than_0.05 <- function(df) {
#   sum(df < 0.05)
# }
# counts <- tapply(data$p_value, data$k, count_p_less_than_0.05)

ggplot(data, aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  facet_wrap(~ group) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.00, 1.00)) +
  labs(x = "p-value", y = "Count", title = "Histogram of p-values for different P_mis")













# library(forestplot)
# library(ggplotify)
# library(patchwork)
# data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
# data <- na.omit(data)
# head(data)
#
# # selected_genes <- c('ENSG00000214226.8', 'ENSG00000116205.12',
# #                     'ENSG00000166002.6', 'ENSG00000219200.11')
# selected_genes <- data$gene[which(eqtlbmaa > 0.95 & PRPP < 0.001)[c(1, 2, 3, 6)]]
# selected_genes
# plots <- list()
# for (i in 1:4) {
#   selected_gene <- data[data$gene == selected_genes[i],]
#
#   plot_matrix <- matrix(c(
#     selected_gene$bhat_Artery_Aorta, selected_gene$bhat_Liver, selected_gene$bhat_Muscle_Skeletal,
#     selected_gene$bhat_Artery_Aorta - 1.96 * selected_gene$se_Artery_Aorta,
#     selected_gene$bhat_Liver - 1.96 * selected_gene$se_Liver,
#     selected_gene$bhat_Muscle_Skeletal - 1.96 * selected_gene$se_Muscle_Skeletal,
#     selected_gene$bhat_Artery_Aorta + 1.96 * selected_gene$se_Artery_Aorta,
#     selected_gene$bhat_Liver + 1.96 * selected_gene$se_Liver,
#     selected_gene$bhat_Muscle_Skeletal + 1.96 * selected_gene$se_Muscle_Skeletal
#   ), nrow=3, byrow=TRUE)
#
#   # Create labels for the plot
#   labels <- c('Artery Aorta', 'Liver', 'Muscle Skeletal')
#
#   p <-
#     forestplot(labeltext=labels,
#                mean=plot_matrix[1,],
#                lower=plot_matrix[2,],
#                upper=plot_matrix[3,],
#                xlab=paste("Effect Size for Gene", selected_genes[i]),
#                lwd.ci = 3, # Make the lines bolder
#                txt_gp = fpTxtGp(
#                  label = gpar(cex = 1.2),      # 调整y轴标签字体大小
#                  ticks = gpar(cex = 1.2),      # 调整x轴刻度标签字体大小
#                  xlab = gpar(cex = 1.4)        # 调整x轴标签字体大小
#                ),
#                boxsize = 0.15, # 调整图形元素大小
#                lineheight = unit(1.5, "cm"), # 增加行高
#                # col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
#                graphwidth = unit(80, "mm"), # 增加图形部分的宽度
#                new_page = F)
#
#   plots[[i]] <- p
# }
#
# p1_grob <- grid2grob(print(plots[[1]]))
# p2_grob <- grid2grob(print(plots[[2]]))
# p3_grob <- grid2grob(print(plots[[3]]))
# p4_grob <- grid2grob(print(plots[[4]]))
#
# p_combined <- wrap_elements(p1_grob) + wrap_elements(p2_grob) +
#   wrap_elements(p3_grob) + wrap_elements(p4_grob)
# p_combined
