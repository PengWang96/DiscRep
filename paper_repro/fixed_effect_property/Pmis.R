rm(list = ls())
library(ggplot2)
# mu1 <- 0
# mu2 <- 1
x <- seq(-10, 10, length=1000)

P_mis <- function(k) {
  integrand <- function(z) {
    1 / (1 + exp(0.5 * (1/k)^2 + z/k)) * dnorm(z)
  }
  result <- integrate(integrand, lower = -Inf, upper = Inf)
  return(result$value)
}

k <- seq(0.1, 4, 0.1); P_mis_values <- sapply(k, P_mis)
plot(k, P_mis_values, type='l', xlab='k', ylab='P_mis', main='P_mis function')
# 
# 
# for (i in seq_along(P_mis_values)) {
#   formatted_P_mis <- sprintf("%.4f", P_mis_values[i])
#   formatted_k <- sprintf("%.4f", k[i])
#   
#   y1 <- dnorm(x, mean = mu1, sd = k[i] * mu2)
#   y2 <- dnorm(x, mean = mu2, sd = k[i] * mu2)
#   
#   df <- data.frame(x = x, y1 = y1, y2 = y2)
#   p <- ggplot(df, aes(x)) +
#     geom_line(aes(y = y1, colour = "Distribution 1"), size = 1) +
#     geom_line(aes(y = y2, colour = "Distribution 2"), size = 1) +
#     geom_ribbon(aes(ymin = pmin(y1, y2), ymax = 0), fill = "purple", alpha = 0.2) +
#     labs(title = paste('Distinguishability Criterion (', 'Pmis: ', formatted_P_mis, ', k:', formatted_k, ')'),
#          x = 'Value', y = 'Density') +
#     theme_minimal() +
#     scale_colour_manual("", 
#                         breaks = c("Distribution 1", "Distribution 2"),
#                         values = c("blue", "red"))
#   
#   print(p)
#   Sys.sleep(1)
# }



rm(list = ls())
library(ggplot2)
library(gridExtra)
mu1 <- 0
mu2 <- 1
x <- seq(-10, 10, length=1000)

P_mis <- function(k) {
  integrand <- function(z) {
    1 / (1 + exp(0.5 * (1/k)^2 + z/k)) * dnorm(z)
  }
  result <- integrate(integrand, lower = -Inf, upper = Inf)
  return(result$value)
}
k_values <- c(0.1, 1, 4)
P_mis_values <- sapply(k_values, P_mis)

plots <- list()

for (i in 1:length(P_mis_values)) {
  formatted_P_mis <- sprintf("%.4f", P_mis_values[i])
  formatted_k <- sprintf("%.4f", k_values[i])
  
  y1 <- dnorm(x, mean = mu1, sd = k_values[i] * mu2)
  y2 <- dnorm(x, mean = mu2, sd = k_values[i] * mu2)
  
  df <- data.frame(x = x, y1 = y1, y2 = y2)
  p <- ggplot(df, aes(x)) +
    geom_line(aes(y = y1, colour = "Distribution 1"), size = 1) +
    geom_line(aes(y = y2, colour = "Distribution 2"), size = 1) +
    geom_ribbon(aes(ymin = pmin(y1, y2), ymax = 0), fill = "purple", alpha = 0.2) +
    labs(title = paste('Pmis: ', formatted_P_mis, ', k:', formatted_k),
         x = 'Value', y = 'Density') +
    theme_minimal() +
    scale_colour_manual("", 
                        breaks = c("Distribution 1", "Distribution 2"),
                        values = c("blue", "red"))
  plots[[i]] <- p
}
grid.arrange(grobs = plots, ncol = 3)













rm(list = ls())
library(ggplot2)
library(gridExtra)
mu1 <- 0
mu2 <- 1
x <- seq(-10, 10, length=1000)

P_mis <- function(k) {
  integrand <- function(z) {
    1 / (1 + exp(0.5 * (1/k)^2 + z/k)) * dnorm(z)
  }
  result <- integrate(integrand, lower = -Inf, upper = Inf)
  return(result$value)
}
k_values <- c(0.1, 1, 4)
P_mis_values <- sapply(k_values, P_mis)

plots <- list()
for (i in 1:length(P_mis_values)) {
  # formatted_P_mis <- sprintf("%.1f", P_mis_values[i])
  # formatted_k <- sprintf("%.4f", k_values[i])
  formatted_P_mis <- c(0, 0.4, 0.5)
  formatted_k <- c(0.1, 1, 4)
  
  y1 <- dnorm(x, mean = mu1, sd = k_values[i] * mu2)
  y2 <- dnorm(x, mean = mu2, sd = k_values[i] * mu2)
  
  df <- data.frame(x = x, y1 = y1, y2 = y2)
  p <- ggplot(df, aes(x)) +
    geom_line(aes(y = y1, colour = "Study 1"), size = 1) +
    geom_line(aes(y = y2, colour = "Study 2"), size = 1) +
    geom_ribbon(aes(ymin = pmin(y1, y2), ymax = 0), fill = "purple", alpha = 0.2) +
    labs(title = paste('Pmis:', formatted_P_mis[i], ', k:', formatted_k[i]),
         x = 'Value', y = 'Density') +
    theme_minimal() +
    scale_colour_manual("", 
                        breaks = c("Study 1", "Study 2"),
                        values = c("blue", "red"))
  if(i == 3) {
    p <- p + theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }
  plots[[i]] <- p
}

do.call(grid.arrange, c(plots, ncol = 3))
