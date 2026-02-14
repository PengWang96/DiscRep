#' Compute Posterior−PRP under distinguishability criterion using Metropolis-Hastings MCMC Algorithm
#'
#' This function computes the posterior probability under the distinguishability criterion
#' using the Metropolis-Hastings algorithm for Markov Chain Monte Carlo (MCMC) simulations.
#' It is designed to test hypotheses of replicability with tolerable heterogeneity levels.
#'
#' @param N Integer. Total number of iterations in the Metropolis-Hastings algorithm.
#' @param r Numeric. Burn-in ratio, representing the proportion of iterations to discard as burn-in (e.g., `0.05`).
#' @param m Integer. Number of replications for test statistics.
#' @param hat_beta Numeric vector. The original or simulated estimated effects.
#' @param hat_sigma_sq Numeric vector. The squared standard errors corresponding to the estimated effects in `hat_beta`.
#' @param test Character. Type of test statistic to use. Options are `"Q"` (heterogeneity test) or `"Egger"` (Egger's regression test).
#' @param side Character. Type of hypothesis test to perform. Options are `"one.side"` (one-tailed) or `"two.side"` (two-tailed).
#' @param k_vec Numeric vector (optional). A precomputed vector of `k` values, representing heterogeneity levels corresponding to a range of misclassification probabilities (`P_mis`) from 0 to 0.05.
#'              If `k_vec` is `NULL`, a default vector is generated using
#'              `pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05); k_vec <- sapply(pvec, inverse_P_mis)`,
#'              and cached for reuse within the current R process.
#'              If `k_vec` has length 1, a fixed heterogeneity level is used.
#'              If `k_vec` has length > 1, `k` is sampled from `k_vec` during MCMC.
#'              Example generation: `pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05); k_vec <- sapply(pvec, inverse_P_mis)`.
#' @param bar_beta_init Numeric scalar (optional). Initial value of `bar_beta` at iteration 1.
#'              If `NULL`, the default initialization is `hat_beta[1]`.
#' @return A list containing the following elements:
#' \item{p_value}{Numeric. The computed p-value for the test.}
#' \item{bar_beta}{Numeric vector. The estimated posterior means for each iteration of MCMC.}
#' \item{acception}{Numeric. The acceptance rate of the Metropolis-Hastings algorithm.}
#' @details
#' This function implements a Metropolis-Hastings algorithm to sample from the posterior
#' distribution of the model parameters under the distinguishability criterion.
#' If `k_vec` is not provided, a default `k_vec` is generated and cached.
#' A warning is issued if any value in `k_vec` is greater than 0.2726814,
#' because only `0 <= k <= 0.2726814` corresponds to `P_mis < 0.05`.
#'
#' @importFrom stats rnorm runif
#' @examples
#' # Generate k_vec using a logarithmic sequence
#' pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
#' k_vec <- sapply(pvec, inverse_P_mis)
#'
#' # Run the Metropolis-Hastings MCMC
#' result <- metropolis_hastings(
#'   N = 10000,
#'   r = 0.05,
#'   m = 2,
#'   hat_beta = c(0.4, 0.6),
#'   hat_sigma_sq = c(0.1, 0.2),
#'   k_vec = k_vec
#' )
#'
#' # Print results
#' print(result$p_value)
#'
#' # Use default cached k_vec
#' result_default <- metropolis_hastings(
#'   N = 10000,
#'   r = 0.05,
#'   m = 2,
#'   hat_beta = c(0.4, 0.6),
#'   hat_sigma_sq = c(0.1, 0.2)
#' )
#'
#' # Set custom bar_beta initialization
#' result_custom_init <- metropolis_hastings(
#'   N = 10000,
#'   r = 0.05,
#'   m = 2,
#'   hat_beta = c(0.4, 0.6),
#'   hat_sigma_sq = c(0.1, 0.2),
#'   bar_beta_init = 0.5
#' )
#' @export
metropolis_hastings <- function(N, r, m, hat_beta, hat_sigma_sq, test = "Q",
                                side = "one.side", k_vec = NULL, bar_beta_init = NULL) {

  if (is.null(k_vec)) {
    k_vec <- get_default_k_vec()
  }

  if (!is.numeric(k_vec) || length(k_vec) == 0 || anyNA(k_vec)) {
    stop("k_vec must be a non-empty numeric vector without NA.")
  }

  if (any(k_vec < 0)) {
    stop("k_vec must be non-negative.")
  }

  if (!is.null(bar_beta_init) &&
      (!is.numeric(bar_beta_init) || length(bar_beta_init) != 1 ||
       is.na(bar_beta_init) || !is.finite(bar_beta_init))) {
    stop("bar_beta_init must be a single finite numeric value.")
  }

  k_pmis_0_05 <- 0.2726814
  if (any(k_vec > k_pmis_0_05)) {
    warning(
      "Use heterogeneity levels > 0.2726814 with caution: only 0 <= k <= 0.2726814 corresponds to P_mis < 0.05.",
      call. = FALSE
    )
  }

  is_fixed_k <- length(k_vec) == 1

  # Initialization
  bar_beta_start <- if (is.null(bar_beta_init)) hat_beta[1] else bar_beta_init
  k_max <- if (is_fixed_k) k_vec[1] else k_pmis_0_05
  # k = 100, 0.5420572, 0.4626489, 0.3349839, 0.2726814 corresponding to P_mis = 0.5, 0.25, 0.2, 0.1, 0.05
  bar_beta <- rep(bar_beta_start, N)
  k <- rep(k_max, N)
  count <- 0
  count2 <- 0
  acception <- 0
  hat_beta_star_matrix <- matrix(nrow = m, ncol = N)  # to store all hat_beta_star values

  # Main loop of the algorithm
  for (i in 2:N) {
    # # Sample bar_beta*
    bar_beta_star <- rnorm(1,
                           mean = mean(hat_beta),
                           sd = sqrt(k_max^2 * mean(hat_beta)^2 + mean(hat_sigma_sq)))


    # Sample k*
    k_star <- if (is_fixed_k) k_max else sample(k_vec, 1)

    # Calculate acceptance ratio a
    a1 <- g(bar_beta_star, hat_beta, hat_sigma_sq, k_star)
    a2 <- q(bar_beta[i-1], hat_beta, hat_sigma_sq, k_max)
    a3 <- g(bar_beta[i-1], hat_beta, hat_sigma_sq, k[i-1])
    a4 <- q(bar_beta_star, hat_beta, hat_sigma_sq, k_max)
    a <- exp(a1 + a2 - a3 - a4)
    a <- min(1, a)

    # Accept or reject
    if (runif(1) < a) {
      bar_beta[i] <- bar_beta_star
      k[i] <- k_star
      acception <- acception + 1
    } else {
      bar_beta[i] <- bar_beta[i-1]
      k[i] <- k[i-1]
    }

    # Calculate test statistic and p-value
    if (i > N * r) {
      # varr <- 1 / (1 / (k[i]^2 * bar_beta[i]^2) + 1 / hat_sigma_sq)
      # meann <- varr * (1 / (k[i]^2 * bar_beta[i]) + hat_beta / hat_sigma_sq)
      if (k[i] != 0) {
        varr <- (k[i]^2 * bar_beta[i]^2 * hat_sigma_sq) / (k[i]^2 * bar_beta[i]^2 + hat_sigma_sq)
        meann <- (hat_sigma_sq * bar_beta[i] + k[i]^2 * bar_beta[i]^2 * hat_beta) / (hat_sigma_sq + k[i]^2 * bar_beta[i]^2)
        # meann <- varr * (k[i]^2 * bar_beta[i] * hat_beta + hat_sigma_sq) / (k[i]^2 * bar_beta[i] * hat_sigma_sq)##############
        beta_star <- rnorm(m,
                           mean = meann,
                           sd = sqrt(varr))
      } else {
        beta_star <- bar_beta[i]
      }


      hat_beta_star <- rnorm(m, mean = beta_star, sd = sqrt(hat_sigma_sq))
      hat_beta_star_matrix[, i] <- hat_beta_star  # store the hat_beta_star values

      if (test == "Q") {
        statistic_star <- T_Q(hat_beta_star, bar_beta[i], k[i], hat_sigma_sq)
        statistic_ori <- T_Q(hat_beta, bar_beta[i], k[i], hat_sigma_sq)
      } else if (test == "Egger") {
        statistic_star <- egger(hat_beta_star, bar_beta[i], k[i], hat_sigma_sq, m)
        statistic_ori <- egger(hat_beta, bar_beta[i], k[i], hat_sigma_sq, m)
      }

      if (test != "Hotelling") {
        if (side == "one.side") {
          if (statistic_star >= statistic_ori) {
            count <- count + 1
          }
        } else {
          if (statistic_star >= statistic_ori) {
            count <- count + 1
          } else {
            count2 <- count2 + 1
          }
        }
      }
    }
  }

  if (test != "Hotelling") {
    if (side == "one.side") {
      p_value <- count / (N - N * r)
    } else {
      p_value <- 2 * min(count/(count + count2), count2/(count + count2))
    }
  } else {
    beta_prime_H <- rowMeans(hat_beta_star_matrix[, (N*r+1):N])  # calculate the mean of hat_beta_star values
    p_value <- hotelling_t_squared(beta_prime_H = beta_prime_H, hat_beta, hat_sigma_sq, df = N-N*r)
  }

  list(p_value = p_value, bar_beta = bar_beta, acception = acception/N)
}




#' Fit Fixed Effect Model using Monte Carlo Algorithm
#'
#' This function performs a fixed effect analysis.
#' It calculates the p-value for the fixed effect model by comparing the test statistic
#' of the observed data with that of the simulated data.
#'
#' @param N Total number of iterations for the Monte Carlo simulation.
#' @param m Number of replications for the test statistic.
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#'
#' @return A list containing the p-value, bar_beta estimates.
#'
#' @examples
#' fixed_effect(10000, 2, c(0.4, 0.6), c(0.1, 0.2))
#'
#' @export
fixed_effect <- function(N, m, hat_beta, hat_sigma_sq) {
  bar_beta <- rep(0, N)
  k <- rep(0, N)
  count <- 0

  for (i in 1:N) {
    varrr <- 1 / sum(1 / hat_sigma_sq)
    meannn <- varrr * sum(hat_beta / hat_sigma_sq)
    bar_beta[i] <- rnorm(1, mean = meannn, sd = sqrt(varrr))

    hat_beta_star <- rnorm(m, mean = bar_beta[i], sd = sqrt(hat_sigma_sq))

    T_Q_star <- T_Q(hat_beta_star, bar_beta[i], k[i], hat_sigma_sq)
    T_Q_ori <- T_Q(hat_beta, bar_beta[i], k[i], hat_sigma_sq)

    if (T_Q_star >= T_Q_ori) {
      count <- count + 1
      }
  }

  p_value <- count / N

  list(p_value = p_value, bar_beta = bar_beta)
}






#' Calculate Posterior-PRP using numerical integration in fixed effect model
#'
#' This function calculates the Posterior-PRP (Posterior Predictive Replication Probability) value
#' using numerical integration for a given set of observed effects and variances. It integrates over
#' the posterior distribution of the true effect sizes `bar_beta` and computes the likelihood of
#' observing `hat_beta_prime` under certain conditions.
#'
#' @param hat_beta A numeric vector representing observed effect sizes.
#' @param hat_sigma_sq A numeric vector representing observed variances. It must have the same length as `hat_beta`.
#' @param tol The maximum tolerance in integration, default 1e-6.
#' @return A list containing:
#'   \item{Integral_Value}{The computed Posterior-PRP value.}
#'   \item{Integration_Error}{The estimated numerical integration error.}
#'
#' @importFrom stats integrate pchisq pf
#' @examples
#' # Example usage
#' hat_beta <- c(1, 1)
#' hat_sigma_sq <- c(0.1, 0.1)
#' calc_PRP_FE_NumInt(hat_beta, hat_sigma_sq)
#'
#' @export
calc_PRP_FE_NumInt <- function(hat_beta, hat_sigma_sq, tol = 1e-6) {
  if (length(hat_sigma_sq) != length(hat_beta)) {
    stop("The lengths of hat_beta and hat_sigma_sq must match.")
  }

  # Number of studies
  m <- length(hat_beta)

  # Calculate weights w_j = 1 / hat_sigma_j^2
  w <- 1 / hat_sigma_sq

  # Calculate posterior mean and variance for bar_beta
  V_bar_beta <- 1 / sum(w)  # Variance
  mu_bar_beta <- V_bar_beta * sum(w * hat_beta)  # Mean

  # Constant c_0
  c_0 <- sum(w * (hat_beta - mu_bar_beta)^2)

  # Define the PDF for the Chi-square distribution with 1 df
  d_chi2_1 <- function(x) {
    ifelse(x >= 0, 1 / sqrt(2 * pi) * exp(-x / 2) / sqrt(x), 0)
  }

  # Define the CDF for the Chi-square distribution with m df
  F_chi2_m <- function(x, m) {
    pchisq(x, df = m)
  }

  # Define the integrand based on your theoretical derivation
  integrand <- function(s) {
    # Chi-squared CDF for m degrees of freedom
    (1 - F_chi2_m(s + c_0, m)) * d_chi2_1(s)
  }
  # Numerical integration using integrate function
  result <- integrate(integrand, lower = 0, upper = Inf, rel.tol = tol)
  # Return the results
  list(Integral_Value = result$value, Integration_Error = result$abs.error)
}






#' Calculate the weighted mean of beta estimates
#'
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @return A numeric value representing the weighted mean of the effect estimates.
#' @examples
#' hat_beta <- c(0.5, 1.0, 1.5)
#' hat_sigma_sq <- c(0.2, 0.3, 0.4)
#' mean_Q(hat_beta, hat_sigma_sq)
#'
#' @export
mean_Q <- function(hat_beta, hat_sigma_sq) {
  w <- 1 / hat_sigma_sq
  return(sum(w * hat_beta) / sum(w))
}

#' Calculate the p-value for Cochran's Q test
#'
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @return A numeric value representing the p-value for Cochran's Q test.
#' @importFrom stats pchisq
#' @examples
#' hat_beta <- c(0.5, 1.0, 1.5)
#' hat_sigma_sq <- c(0.2, 0.3, 0.4)
#' frequency_pvalue(hat_beta, hat_sigma_sq)
#' @export
frequency_pvalue <- function(hat_beta, hat_sigma_sq) {
  statistic <- T_Q(hat_beta, mean_Q(hat_beta, hat_sigma_sq), 0, hat_sigma_sq)
  p_value <- 1 - pchisq(statistic, df = length(hat_beta) - 1)
  return(p_value)
}




#' Probability of Misclassification Function
#'
#' This function calculates the probability of misclassification for a given value of k.
#' @param k A numeric value representing the level of heterogeneity.
#' @return Numeric value of the probability of misclassification.
#' @importFrom stats dnorm
#' @examples
#' P_mis(1)
#' @export
P_mis <- function(k) {
  if (k == 0) {
    return(0)
  }
  integrand <- function(z) {
    1 / (1 + exp(0.5 * (1/k)^2 + z/k)) * dnorm(z)
  }
  result <- integrate(integrand, lower = -Inf, upper = Inf)
  return(result$value)
}




#' Inverse Probability of Misclassification Function
#'
#' This function approximates the inverse function value of \eqn{P_{\text{mis}}(k)} using the bisection method.
#' @param P_mis_val The value of misclassification probability for which the inverse is sought.
#' @param lower The lower bound for the bisection method.
#' @param upper The upper bound for the bisection method.
#' @param tol The tolerance level for the bisection method.
#' @param max_iter The maximum number of iterations for the bisection method.
#' @return Numeric value of the approximate inverse function value of \eqn{P_{\text{mis}}}.
#' @examples
#' inverse_P_mis(0.05)
#' @export
inverse_P_mis <- function(P_mis_val, lower = 1e-6, upper = 100,
                          tol = 1e-10, max_iter = 10^15) {
  root_func <- function(k) P_mis(k) - P_mis_val
  iter <- 0
  while (upper - lower > tol && iter < max_iter) {
    mid <- (lower + upper) / 2
    if (root_func(mid) * root_func(lower) <= 0) {
      upper <- mid
    } else {
      lower <- mid
    }
    iter <- iter + 1
  }
  return((lower + upper) / 2)
}

# Cache default k-grid once per R process to avoid repeated inverse_P_mis calls.
.disc_rep_cache <- new.env(parent = emptyenv())

get_default_k_vec <- function() {
  k_vec_default <- get0("k_vec_default", envir = .disc_rep_cache, inherits = FALSE)
  if (is.null(k_vec_default)) {
    pvec <- c(10^seq(-10, log10(0.05), 0.01), 0.05)
    k_vec_default <- sapply(pvec, inverse_P_mis)
    assign("k_vec_default", k_vec_default, envir = .disc_rep_cache)
  }
  k_vec_default
}








#' Function g by log transformation
#'
#' This function calculates the log-transformed value of g.
#' @param bar_beta The estimated true underlying effect.
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @param k The level of heterogeneity.
#' @return Numeric value of the log-transformed function g.
#' @keywords internal
g <- function(bar_beta, hat_beta, hat_sigma_sq, k) {
  m <- length(hat_beta)  # Number of estimated effects
  sum_log_terms <- 0  # Initialize sum of log terms

  for (j in 1:m) {
    sigma_j_sq <- hat_sigma_sq[j]
    beta_j <- hat_beta[j]

    # Calculate the variance term
    variance_term <- k^2 * bar_beta^2 + sigma_j_sq

    # Log of the first part of the product (normalization factor)
    log_normalization <- -0.5 * log(variance_term)

    # Log of the second part (exponent term)
    log_exponent <- -0.5 * ((beta_j - bar_beta)^2 / variance_term)

    # Add these log terms
    sum_log_terms <- sum_log_terms + log_normalization + log_exponent
  }
  return(sum_log_terms)
}





#' Proposed Function q by log transformation
#'
#' This function calculates logarithm of the proposed density used in Metropolis-Hastings (MH).
#' @param bar_beta_star The simulated value of bar_beta.
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @param k The level of heterogeneity.
#' @return Logarithm of the density of the normal distribution for bar_beta_star.
#' @importFrom stats dnorm
#' @keywords internal
q <- function(bar_beta_star, hat_beta, hat_sigma_sq, k) {
  k_max <- k
  res <- dnorm(bar_beta_star,
               mean = mean(hat_beta),
               sd = sqrt(k_max^2 * mean(hat_beta)^2 + mean(hat_sigma_sq)), log = T)
  return(res)
}


# q <- function(bar_beta_star, hat_beta, hat_sigma_sq, k) {
#   varrr <- 1 / sum(1 / hat_sigma_sq)
#   meannn <- varrr * sum(hat_beta / hat_sigma_sq)
#   res <- dnorm(bar_beta_star,
#                mean = meannn,
#                sd = sqrt(varrr), log = T)
#   return(res)
# }







#' Modified Cochran’s Q statistic T_Q
#'
#' This function calculates the modified Cochran’s Q statistic.
#' @param hat_beta The original or the simulated estimated effects.
#' @param bar_beta The estimated true underlying effect.
#' @param k The level of heterogeneity.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @return Numeric value of the test statistic T_Q.
#' @examples
#' T_Q(c(0.4, 0.6), 0.5, 1, c(0.1, 0.2))
#' @export
T_Q <- function(hat_beta, bar_beta, k, hat_sigma_sq) {
  w <- 1 / (hat_sigma_sq + k^2 * bar_beta^2)
  sum(w * (hat_beta - bar_beta)^2)
}




#' Modified Egger test statistics
#'
#' This function provides the calculation for Egger test quantities.
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @param bar_beta The estimated true underlying effect.
#' @param k The level of heterogeneity.
#' @param m The number of replications.
#' @return The Egger test statistic value.
#' @examples
#' egger(c(0.4, 0.6), 0.5, 1, c(0.1, 0.2), 2)
#' @export
egger <- function(hat_beta, bar_beta, k, hat_sigma_sq, m){
  y <- hat_beta/sqrt(hat_sigma_sq + k^2*bar_beta^2)
  x <- 1/sqrt(hat_sigma_sq + k^2*bar_beta^2)

  Sxx <- sum( (x-mean(x))*(x-mean(x)) )
  Sxy <- sum( (x-mean(x))*(y-mean(y)) )
  Syy <- sum( (y-mean(y))*(y-mean(y)) )

  b1 <- Sxy/Sxx
  b0 <- mean(y)- b1*mean(x)

  s2 <- (Syy - b1^2*Sxx)/ (m - 2)
  vb0 <- s2 * (1 / m + mean(x)^2 / Sxx)

  egger_test <- b0^2/vb0
  return(egger_test)
}




#' Hotelling T-squared Statistic Function
#'
#' This function calculates the Hotelling T-squared statistic for multivariate hypothesis testing.
#' @param beta_prime_H The mean of estimated true underlying effect.
#' @param hat_beta The original or the simulated estimated effects.
#' @param hat_sigma_sq The squared standard errors of the estimated effects.
#' @param df Degrees of freedom for the test.
#' @return p-value for the Hotelling T-squared test.
#' @importFrom stats pf
#' @keywords internal
hotelling_t_squared <- function(beta_prime_H, hat_beta, hat_sigma_sq, df) {
  S <- diag(hat_sigma_sq)
  m <- length(hat_beta)
  T_H <- df * t(beta_prime_H - hat_beta) %*% solve(S) %*% (beta_prime_H - hat_beta)
  F_value <- (df - m) / ((df - 1) * m) * T_H
  p_value <- 1 - pf(F_value, m, df - m)
  return(p_value)
}



#' Combine Results Function
#'
#' This function combines the results from multiple Metropolis-Hastings simulations.
#' @param ... Variable number of Metropolis-Hastings result lists.
#' @return A list containing combined p-values, bar_beta estimates, and average acceptance rate.
#' @examples
#' result1 <- list()
#' result2 <- list()
#' result3 <- list()
#' combine_results(result1, result2, result3)
#' @export
combine_results <- function(...) {
  args <- list(...)

  combined_p_value <- c()
  combined_bar_beta <- NULL
  combined_acception <- 0

  for (arg in args) {
    combined_p_value <- c(combined_p_value, arg$p_value)
    combined_acception <- combined_acception + arg$acception

    if (is.null(combined_bar_beta)) {
      combined_bar_beta <- arg$bar_beta
    } else {
      combined_bar_beta <- rbind(combined_bar_beta, arg$bar_beta)
    }
  }

  total_acception <- combined_acception / length(args)

  list(p_value = combined_p_value,
       bar_beta = combined_bar_beta,
       acception = total_acception)
}



#' Shuffle Function
#'
#' This function randomly shuffles elements in a vector.
#' @param xv Input vector to be shuffled.
#' @param rep Number of times to perform the shuffle.
#' @return Shuffled vector.
#' @examples
#' shuffle(c(1, 2, 3, 4, 5), 2)
#' @export
shuffle <- function(xv, rep){
  yv = xv
  for (i in 1:rep){
    index = sample(1:length(xv),2, replace=F)
    temp = yv[index[1]]
    yv[index[1]] = yv[index[2]]
    yv[index[2]] = temp
  }
  return(yv)
}




#' Simulate Batch Contamination Function
#'
#' This function simulates batch contamination in a dataset by incorporating batch effects and heterogeneity in group data.
#' @param bb_sd Standard deviation for batch effects. This determines the variability of the batch effects in the contaminated groups.
#' @param k Level of heterogeneity, affecting the variability of the true underlying effect across groups.
#' @param bbar True underlying effect. This is the baseline effect size which is perturbed by heterogeneity.
#' @param m Number of groups to simulate.
#' @param bc_grp_num Number of groups that are batch-contaminated.
#' @param sample_size Number of samples per group.
#' @param Gv Matrix of group indicators; each row corresponds to a group and each column to a sample.
#' @param Batch Matrix of batch effects; structured similarly to the \code{Gv} matrix, with each row corresponding to a group.
#' @param rnse Standard deviation of the random noise added to the simulation.
#' @return A vector containing the simulated data, specifically the estimated coefficients and their standard errors from linear models applied to each group.
#' @importFrom stats rnorm lm
#' @examples
#' rnse <- 1
#' sample_size <- 5
#' m <- 1
#' bc_grp_num <- 1 # number of batch-contaminated groups
#' Gv <- t(sapply(1:m, function(x) rbinom(sample_size, 1, 0.4)))
#' Batch <- t(apply(Gv, 1, function(x) shuffle(x, rep=floor(sample_size*0.2))))
#' sim_batch(0.1, 1, 0.5, 1, bc_grp_num, sample_size, Gv, Batch, rnse)
#' @export
sim_batch <- function(bb_sd, k, bbar, m, bc_grp_num, sample_size, Gv, Batch, rnse) {
  grp_num <- m

  btv <- bbar + rnorm(grp_num, sd = sqrt(k^2 * bbar^2))

  bbv <- rep(0, grp_num)
  bbv[1:bc_grp_num] <- rnorm(bc_grp_num, sd = bb_sd)

  rst <- c()
  for (i in 1:grp_num) {
    gv <- Gv[i,]
    batch <- Batch[i,]
    y <- btv[i] * gv + bbv[i] * batch + rnorm(sample_size, sd = rnse)
    dat <- data.frame(y = y, gv = gv)
    fit <- summary(lm(y ~ gv, data = dat))
    rst <- c(rst, fit$coef[2, 1], fit$coef[2, 2])
  }

  return(rst)
}





#' Censoring Function
#'
#' This function applies a censoring mechanism based on a power law.
#' @param p Numeric vector of probabilities to be censored.
#' @param c Censoring strength
#' @return Censored probabilities.
#' @examples
#' wipi(0.1, 1)
#' @export
wipi<-function(p, c){
  return(exp(-c*p^(1.5)))
}
