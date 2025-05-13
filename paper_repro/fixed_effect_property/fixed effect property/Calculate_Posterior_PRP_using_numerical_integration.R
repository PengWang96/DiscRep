#' Calculate Posterior-PRP using numerical integration
#'
#' This function calculates the Posterior-PRP (Posterior Predictive Replication Probability) value
#' using numerical integration for a given set of observed effects and variances. It integrates over
#' the posterior distribution of the true effect sizes `bar_beta` and computes the likelihood of
#' observing `hat_beta_prime` under certain conditions.
#'
#' @param hat_beta A numeric vector representing observed effect sizes.
#' @param hat_sigma_sq A numeric vector representing observed variances. It must have the same length as `hat_beta`.
#' @return A list containing:
#'   \item{Integral_Value}{The computed Posterior-PRP value.}
#'   \item{Integration_Error}{The estimated numerical integration error.}
#' @importFrom cubature hcubature
#' @examples
#' # Example usage
#' if (!requireNamespace("cubature", quietly = TRUE)) {
#' install.packages("cubature")
#' }
#' library(cubature)
#' hat_beta <- c(1, 1)
#' hat_sigma_sq <- c(0.1, 0.1)
#' calculate_posterior_prp(hat_beta, hat_sigma_sq)
#' 
#' @export
calculate_posterior_prp <- function(hat_beta, hat_sigma_sq) {
  if (length(hat_sigma_sq) != length(hat_beta)) {
    stop("The lengths of hat_beta and hat_sigma_sq must match.")
  }
  
  # Inner function to calculate T_Q
  T_Q <- function(hat_beta, bar_beta, hat_sigma_sq) {
    w <- 1 / hat_sigma_sq
    sum(w * (hat_beta - bar_beta)^2)
  }
  
  # Posterior density for bar_beta
  posterior_bar_beta <- function(bar_beta, hat_beta, hat_sigma_sq) {
    varrr <- 1 / sum(1 / hat_sigma_sq)
    meannn <- varrr * sum(hat_beta / hat_sigma_sq)
    dnorm(bar_beta, mean = meannn, sd = sqrt(varrr))
  }
  
  # Likelihood for hat_beta_prime given bar_beta
  likelihood_hat_beta_prime <- function(hat_beta_prime, bar_beta, hat_sigma_sq) {
    dnorm(hat_beta_prime, mean = bar_beta, sd = sqrt(hat_sigma_sq))
  }
  
  # Integrand function
  integrand <- function(bar_beta, hat_beta_prime, hat_beta, hat_sigma_sq) {
    T_prime <- T_Q(hat_beta_prime, bar_beta, hat_sigma_sq)
    T_observed <- T_Q(hat_beta, bar_beta, hat_sigma_sq)
    ifelse(T_prime >= T_observed,
           prod(likelihood_hat_beta_prime(hat_beta_prime, bar_beta, hat_sigma_sq)) *
             posterior_bar_beta(bar_beta, hat_beta, hat_sigma_sq),
           0)
  }
  
  # Define the integrand for cubature
  integrand_cubature_general <- function(x) {
    bar_beta <- x[1]
    hat_beta_prime <- x[-1]
    integrand(bar_beta, hat_beta_prime, hat_beta, hat_sigma_sq)
  }
  
  # Numerical integration using hcubature
  lower_limit <- rep(-Inf, length(hat_beta) + 1)
  upper_limit <- rep(Inf, length(hat_beta) + 1)
  result <- hcubature(integrand_cubature_general, lowerLimit = lower_limit, upperLimit = upper_limit, tol = 1e-2)
  
  
  # Return results
  list(Integral_Value = result$integral, Integration_Error = result$error)
}
