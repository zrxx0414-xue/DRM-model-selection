#generate data
simulate_drm_selection <- function(
    B, correct_model, n_train = 200, n_test = 200, beta_1 = 1, delta = 0,
    model_1_formula = y ~ x1, model_2_formula = y ~ x1 + x2 + x3,
    sigma_error = 1, x1_generator = unit_variance_generators$normal,
    x2_generator = unit_variance_generators$normal, x3_generator = unit_variance_generators$normal,
    data_error_generator = unit_variance_generators$normal,
    test_seed = 999, train_seed = 123) {
  
  stopifnot(correct_model %in% c("model_1", "model_2"))
  generate_predictors <- function(n) {
    x1 <- x1_generator(n)
    x2 <- x2_generator(n)
    x3 <- x3_generator(n)
    data.frame(x1 = x1, x2 = x2, x3 = x3)
  }
  add_response <- function(x, with_noise) {
    signal <- 1 + beta_1 * x$x1 + delta * x$x2 + delta * x$x3
    x$y <- signal + if (with_noise) sigma_error * data_error_generator(nrow(x)) else 0  #add noise in train data
    x
  }
  
  set.seed(test_seed)
  test_data <- add_response(generate_predictors(n_test), with_noise = FALSE)  #generate test data without noise
  set.seed(train_seed)
  iteration_results <- vector("list", B)
  coefficient_results <- vector("list", B)
  
# reapeat the selection process 500 times
  for (b in seq_len(B)) {
    train_x <- generate_predictors(n_train)
    train_without <- add_response(train_x, with_noise = FALSE)
    train_with <- add_response(train_x, with_noise = TRUE)
    y0 <- test_data$y
    
    failed_scores <- function() data.frame(
      metric = drm_metrics, model_1 = NA_real_, model_2 = NA_real_,
      delta_score = NA_real_, selected_model = NA_character_
    )
    analyse_training_condition <- function(train_data, condition) {
      model_1 <- lm(model_1_formula, data = train_data)
      model_2 <- lm(model_2_formula, data = train_data)
      mu_1 <- predict(model_1, newdata = test_data)
      mu_2 <- predict(model_2, newdata = test_data)
      drm <- tryCatch(drm_scores(y0, mu_1, mu_2), error = function(e) failed_scores())
      aic <- data.frame(metric = "AIC", model_1 = AIC(model_1), model_2 = AIC(model_2)) |>
        mutate(delta_score = model_2 - model_1, selected_model = ifelse(model_1 < model_2, "model_1", "model_2"))
      scores <- bind_rows(drm, aic) |> mutate(prediction_type = condition)
      coefficients <- bind_rows(
        data.frame(model = "M1", term = names(coef(model_1)), estimate = unname(coef(model_1))),
        data.frame(model = "M2", term = names(coef(model_2)), estimate = unname(coef(model_2)))
      ) |> mutate(training_condition = condition)
      list(scores = scores, coefficients = coefficients)
    }
    without <- analyse_training_condition(train_without, "Without noise")
    with <- analyse_training_condition(train_with, "With noise")
    iteration_results[[b]] <- bind_rows(without$scores, with$scores) |>
      mutate(iteration = b, correct_model = correct_model, correct = selected_model == correct_model)
    coefficient_results[[b]] <- bind_rows(without$coefficients, with$coefficients) |>
      mutate(iteration = b) |>
      mutate(term = recode(term, `(Intercept)` = "Intercept", x1 = "X1", x2 = "X2", x3 = "X3"))
  }
  
  all_results <- bind_rows(iteration_results)
  coefficients <- bind_rows(coefficient_results) |>
    mutate(term = factor(term, levels = c("Intercept", "X1", "X2", "X3")), model = factor(model, levels = c("M1", "M2")))
  summary <- all_results |>
    group_by(prediction_type, metric) |>
    summarise(

      valid_runs = sum(!is.na(correct)), failed_runs = B - valid_runs,
      failure_rate = failed_runs / B, correct_times = sum(correct, na.rm = TRUE),
      accuracy = ifelse(valid_runs > 0, correct_times / valid_runs, NA_real_),
      mc_se = ifelse(valid_runs > 0, sqrt(accuracy * (1 - accuracy) / valid_runs), NA_real_),
      ci_lower = pmax(0, accuracy - 1.96 * mc_se), ci_upper = pmin(1, accuracy + 1.96 * mc_se),
      mean_delta = mean(delta_score, na.rm = TRUE), median_delta = median(delta_score, na.rm = TRUE),
      delta_lower = quantile(delta_score, 0.025, na.rm = TRUE, names = FALSE),
      delta_upper = quantile(delta_score, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) |>
    mutate(
      metric = factor(metric, levels = metric_order),
      prediction_type = factor(prediction_type, levels = prediction_type_order)
    ) |>
    arrange(prediction_type, metric)
  
  list(summary = summary, all_results = all_results, coefficients = coefficients, test_data = test_data)
}