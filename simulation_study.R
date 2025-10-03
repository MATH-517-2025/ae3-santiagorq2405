# ============================================================================
# MATH-517 Assignment 3 
# Author: Santiago Rivadeneira Quintero
# Date: 03/10/2025
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(knitr)
})

set.seed(2024)

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

# True regression function
m_true <- function(x) {
  sin((x/3 + 0.1)^(-1))
}

# Quartic (biweight) kernel
K_quartic <- function(u) {
  ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0)
}

# Generate data
generate_data <- function(n, alpha = 2, beta = 2, sigma2 = 0.15) {
  if (n < 10) stop("Sample size must be at least 10")
  if (alpha <= 0 || beta <= 0) stop("Beta parameters must be positive")
  if (sigma2 <= 0) stop("Variance must be positive")
  
  X <- rbeta(n, alpha, beta)
  epsilon <- rnorm(n, 0, sqrt(sigma2))
  Y <- m_true(X) + epsilon
  
  data.frame(X = X, Y = Y)
}

# Block method to estimate theta_22 and sigma^2
block_method <- function(data, N) {
  n <- nrow(data)
  
  if (N < 1) stop("N must be at least 1")
  if (n - 5*N <= 0) stop("N too large: need n - 5*N > 0")
  
  data <- data[order(data$X), ]
  
  if (N == 1) {
    data$block <- rep(1, n)
  } else {
    data$block <- cut(rank(data$X, ties.method = "first"), 
                      breaks = N, labels = FALSE, include.lowest = TRUE)
  }
  
  theta_22_sum <- 0
  sigma2_sum <- 0
  
  for (j in 1:N) {
    block_data <- data[data$block == j, ]
    
    if (nrow(block_data) >= 6) {
      model <- lm(Y ~ poly(X, 4, raw = TRUE), data = block_data)
      beta_hat <- coef(model)
      
      # Vectorized operations
      X_block <- block_data$X
      Y_hat <- fitted(model)
      m_j_double_prime <- 2*beta_hat[3] + 6*beta_hat[4]*X_block + 12*beta_hat[5]*X_block^2
      
      theta_22_sum <- theta_22_sum + sum(m_j_double_prime^2)
      sigma2_sum <- sigma2_sum + sum((block_data$Y - Y_hat)^2)
    }
  }
  
  list(
    theta_22 = theta_22_sum / n,
    sigma2 = sigma2_sum / (n - 5*N),
    RSS = sigma2_sum
  )
}

# Compute h_AMISE
compute_h_AMISE <- function(n, sigma2_hat, theta_22_hat, support_length = 1) {
  if (theta_22_hat <= 0) theta_22_hat <- abs(theta_22_hat)
  if (sigma2_hat <= 0) sigma2_hat <- abs(sigma2_hat)
  
  n^(-1/5) * (35 * sigma2_hat * support_length / theta_22_hat)^(1/5)
}

# Mallow's Cp for N selection (Ruppert et al. 1995)
mallows_cp <- function(data, N_values = NULL) {
  n <- nrow(data)
  N_max <- max(min(floor(n / 20), 5), 1)  # CORRECT formula
  
  if (n - 5*N_max <= 0) N_max <- max(floor((n - 1) / 5), 1)
  if (is.null(N_values)) N_values <- 1:N_max
  
  RSS_values <- numeric(length(N_values))
  
  for (i in seq_along(N_values)) {
    N <- N_values[i]
    if (N >= 1 && N <= N_max && n - 5*N > 0) {
      tryCatch({
        RSS_values[i] <- block_method(data, N)$RSS
      }, error = function(e) {
        RSS_values[i] <- NA
      })
    } else {
      RSS_values[i] <- NA
    }
  }
  
  RSS_N_max <- tryCatch(block_method(data, N_max)$RSS, error = function(e) NA)
  
  if (is.na(RSS_N_max) || RSS_N_max <= 0) {
    return(data.frame(N = N_values, Cp = NA, RSS = RSS_values))
  }
  
  Cp_values <- RSS_values / (RSS_N_max / (n - 5*N_max)) - (n - 10*N_values)
  data.frame(N = N_values, Cp = Cp_values, RSS = RSS_values)
}

# Find optimal N
find_optimal_N <- function(data) {
  n <- nrow(data)
  N_max <- max(min(floor(n / 20), 5), 1)  # CORRECT formula
  
  if (n - 5*N_max <= 0) N_max <- max(floor((n - 1) / 5), 1)
  
  cp_results <- mallows_cp(data, 1:N_max)
  cp_results <- cp_results[!is.na(cp_results$Cp), ]
  
  if (nrow(cp_results) > 0) {
    return(cp_results$N[which.min(cp_results$Cp)])
  } else {
    return(1)
  }
}

# ============================================================================
# EXAMPLE: How to run the simulation
# ============================================================================

cat("\nTo run the complete simulation, uncomment the code below:\n\n")

# # Example usage:
# data <- generate_data(n = 500, alpha = 2, beta = 2)
# N_opt <- find_optimal_N(data)
# est <- block_method(data, N_opt)
# h <- compute_h_AMISE(500, est$sigma2, est$theta_22)
# cat("Optimal bandwidth:", h, "\n")

