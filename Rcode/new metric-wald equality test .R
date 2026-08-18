#equality test (WALD)
equality_wald_test <- function(y0, y1, y2) {
  all_y <- c(y0, y1, y2)
  #I fit standardized value in to avoid extreme value.
  centre <- mean(all_y)
  scale_value <- sd(all_y)
  basis_interaction <- function(y) {
    z <- (y - centre) / scale_value
    c(z = z, z2 = z^2, z3 = z^3, z4 = z^4)
  }
  fit <- drmdel(x = all_y, n_samples = c(length(y0), length(y1), length(y2)), basis_func = basis_interaction, control = list(maxit = 2000))
  covariance <- as.matrix(meleCov(fit))
  d <- 4L
  beta_1_index <- 2:(d + 1L)
  beta_2_index <- (d + 3L):(2L * d + 2L)
  difference <- fit$mele[beta_1_index] - fit$mele[beta_2_index]
  
  contrast <- matrix(0, nrow = d, ncol = 2L * (d + 1L))
  for (j in seq_len(d)) {
    contrast[j, beta_1_index[j]] <- 1
    contrast[j, beta_2_index[j]] <- -1
  }
  
  difference_covariance <- contrast %*% covariance %*% t(contrast)
  statistic <- length(all_y) * drop(t(difference) %*% solve(difference_covariance, difference))
  data.frame(wald = statistic, p_value = pchisq(statistic, df = d, lower.tail = FALSE))
}