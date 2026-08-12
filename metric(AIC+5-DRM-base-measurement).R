library(drmdel)

drm_compare <- function(y0, y1, y2, basis_function, grid_length = 200,model_fit_1 = NULL, model_fit_2 = NULL) {
  all_y <- c(y0, y1, y2)
  first_basis <- basis_function(all_y[1])
  d <- length(first_basis)
  
  basis_matrix <- function(y) {
    values <- vapply(y, function(z) as.numeric(basis_function(z)), numeric(d))
    t(matrix(values, nrow = d, ncol = length(y)))
  }
  
  fit <- drmdel(
    x = all_y,
    n_samples = c(length(y0), length(y1), length(y2)),
    basis_func = basis_function,
    control = list(maxit = 1000)
  )
  
  theta <- matrix(fit$mele, nrow = 2, byrow = TRUE)
  
  log_ratio <- function(y, model) {
    Q <- basis_matrix(y)
    as.numeric(theta[model, 1] + Q %*% theta[model, -1])
  }
  
  # Log-ratio distance
  logratio_1 <- mean(log_ratio(y0, 1)^2)
  logratio_2 <- mean(log_ratio(y0, 2)^2)
  
  # Jeffreys divergence
  jeffreys_1 <- mean(log_ratio(y1, 1)) - mean(log_ratio(y0, 1))
  jeffreys_2 <- mean(log_ratio(y2, 2)) - mean(log_ratio(y0, 2))
  
  # CDF distance
  grid <- seq(min(all_y), max(all_y), length.out = grid_length)
  cdf_fit <- cdfDRM(k = 0:2, x = grid, drmfit = fit)
  
  F0 <- cdf_fit[[1]]$cdf_est
  F1 <- cdf_fit[[2]]$cdf_est
  F2 <- cdf_fit[[3]]$cdf_est
  
  cdf_1 <- mean((F1 - F0)^2)
  cdf_2 <- mean((F2 - F0)^2)
  
  # Joint Wald statistics
  covariance <- as.matrix(meleCov(fit))
  n_total <- length(all_y)
  
  beta1_index <- 2:(d + 1)
  beta2_index <- (d + 3):(2 * d + 2)
  
  beta1 <- fit$mele[beta1_index]
  beta2 <- fit$mele[beta2_index]
  
  wald_1 <- tryCatch(
    n_total * drop(t(beta1) %*% solve(
      covariance[beta1_index, beta1_index, drop = FALSE], beta1
    )),
    error = function(e) NA_real_
  )
  
  wald_2 <- tryCatch(
    n_total * drop(t(beta2) %*% solve(
      covariance[beta2_index, beta2_index, drop = FALSE], beta2
    )),
    error = function(e) NA_real_
  )
  
  # Joint DELR: H0 beta1 = 0
  null_fit_1 <- tryCatch(
    drmdel(
      x = all_y,
      n_samples = c(length(y0), length(y1), length(y2)),
      basis_func = basis_function,
      g_null = function(gamma) c(rep(0, d), gamma),
      g_null_jac = function(gamma) rbind(matrix(0, d, d), diag(d)),
      par_dim_null = d,
      control = list(maxit = 1000)
    ),
    error = function(e) NULL
  )
  
  # Joint DELR: H0 beta2 = 0
  null_fit_2 <- tryCatch(
    drmdel(
      x = all_y,
      n_samples = c(length(y0), length(y1), length(y2)),
      basis_func = basis_function,
      g_null = function(gamma) c(gamma, rep(0, d)),
      g_null_jac = function(gamma) rbind(diag(d), matrix(0, d, d)),
      par_dim_null = d,
      control = list(maxit = 1000)
    ),
    error = function(e) NULL
  )
  
  delr_1 <- if (is.null(null_fit_1)) NA_real_ else max(0, as.numeric(null_fit_1$delr))
  delr_2 <- if (is.null(null_fit_2)) NA_real_ else max(0, as.numeric(null_fit_2$delr))
  
  # AIC of candidate regression models
  aic_1 <- if (is.null(model_fit_1)) NA_real_ else tryCatch(AIC(model_fit_1), error = function(e) NA_real_)
  aic_2 <- if (is.null(model_fit_2)) NA_real_ else tryCatch(AIC(model_fit_2), error = function(e) NA_real_)
  
  result <- data.frame(
    metric = c("log_ratio", "Jeffreys", "CDF", "Joint_Wald", "Joint_DELR", "AIC"),
    model_1 = c(logratio_1, jeffreys_1, cdf_1, wald_1, delr_1, aic_1),
    model_2 = c(logratio_2, jeffreys_2, cdf_2, wald_2, delr_2, aic_2),
    stringsAsFactors = FALSE
  )
  
  result$selected_model <- mapply(function(value_1, value_2) {
    if (!is.finite(value_1) || !is.finite(value_2)) return(NA_character_)
    if (isTRUE(all.equal(value_1, value_2))) return("tie")
    if (value_1 < value_2) "model_1" else "model_2"
  }, result$model_1, result$model_2, USE.NAMES = FALSE)
  
  result
}