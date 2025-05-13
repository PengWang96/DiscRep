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


library(forestplot)
library(ggplotify)
library(patchwork)
data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
data <- na.omit(data)
head(data)
selected_genes <- data$gene[which(eqtlbmaa > 0.7 & PRPP < 0.05)] # [c(5, 8, 9, 11, 12)]




dataaa <- data[which(data$gene %in% selected_genes), ]
a <- which(rowSums(dataaa[, c(2, 4, 6)] > 0) == 3)
b <- which(rowSums(dataaa[, c(2, 4, 6)] < 0) == 3)
c <- dataaa[c(a,b), ]
c$gene
# ENSG00000112877.7 ENSG00000108784.9  ENSG00000162437.14  ENSG00000119242.8

for (ii in 1:length(c$gene)) {
  selected_genes <- "ENSG00000112877.7" # c$gene[ii]
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

  # par(mar = c(1, 1, 1, 1))
  p <- forestplot(labeltext=labels,
                  mean=plot_matrix[1,],
                  lower=plot_matrix[2,],
                  upper=plot_matrix[3,],
                  xlab=paste("Effect Size for Gene", selected_genes[gene_index]),
                  lwd.ci = 4, # Make the lines bolder
                  txt_gp = fpTxtGp(
                    label = gpar(cex = 1.6),    # 调整y轴标签字体大小
                    ticks = gpar(cex = 1.2),    # 调整x轴刻度标签字体大小
                    xlab = gpar(cex = 1.6)      # 调整x轴标签字体大小
                  ),
                  boxsize = 0.12,               # 调整图形元素大小
                  lineheight = unit(6, "cm"), # 增加行高
                  graphwidth = unit(150, "mm"),  # 增加图形部分的宽度
                  # col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
                  new_page = T)
  print(p)
}


