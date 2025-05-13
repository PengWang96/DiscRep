summ <- rowSums(data[, 2:8]) %>% as.data.frame()
hist(summ$.)



library(vroom)
data <- vroom("./data/posterior_pvalue.txt.gz",
              delim = "\t",
              skip = 1,
              col_types = cols(
                .default = col_double(),
                gene = col_character(),
                snp = col_character()
              ),
              num_threads = 32)
print(head(data))



library(vroom)
data <- vroom("./data/marginal_pvalue.txt.gz",
              skip = 0,
              col_select = c(1, 2, 3, 10),
              col_types = cols(
                .default = col_double(),
                gene = col_character(),
                snp = col_character()
              ),
              num_threads = 32)
pval <- 1 - data[, 3] + data[, 4]
hist(pval[[1]])
sum(pval[[1]] < 0.05)


data <- vroom("./data/output_prefix_avg_bfs.txt.gz",
              skip = 0,
              delim = "\t",
              n_max = 10000,
              col_types = cols(
                .default = col_double(),
                gene = col_character(),
                snp = col_character()
              ),
              num_threads = 32)


data <- vroom("./data/out_eqtlbma_avg_bfs_genes.txt.gz",
              skip = 0,
              delim = "\t",
              n_max = 10000,
              num_threads = 32)
library(dplyr)
dataa <- data %>%
  filter(gene == "ENSG00000001497.16") %>%
  mutate(sum(snp.post.the))

data <- vroom("./data/out_eqtlbma_avg_bfs_genes.txt.gz",
              skip = 0,
              delim = "\t",
              # n_max = 1000,
              num_threads = 32)



a <- result$p_values
b <- result$p_values[1:1000]
plot(a, b)



rm(list = ls())
setwd("D:/R/Replicability/Peng/real data/eqtlbma/different_tissuses_combination/three_tissues")
library(vroom)
library(dplyr)
data <- vroom("./data/output_prefix_avg_bfs.txt.gz",
              skip = 1,
              delim = "\t",
              n_max = 500000,
              # n_max = 50,
              col_types = cols(
                .default = col_double(),
                gene = col_character(),
                snp = col_character()
              ),
              num_threads = 32)
head(data)

target_genes <- c("ENSG00000003402.19", "ENSG00000001497.16", "ENSG00000004864.13")

filtered_data <- data %>% filter(gene %in% target_genes)
filtered_data <- filtered_data %>%
  mutate(
    sum_post_config = rowSums(select(., post.config.1, post.config.2, post.config.3, 
                                     `post.config.1-2`, `post.config.1-3`, `post.config.2-3`, 
                                     `post.config.1-2-3`), na.rm = TRUE)
  )

result <- filtered_data %>%
  group_by(gene) %>%
  summarize(sum_snp_post_the = sum(snp.post.the, na.rm = TRUE))
sum(filtered_data$sum_post_config < 0.99)
hist(filtered_data$sum_post_config)
