get_ci_params = function(m){
    m_sum = summary(m)
    coef = as.data.frame(m_sum$coefficients)
    coef$low = coef$Estimate - 1.96 * coef$`Std. Error`
    coef$high = coef$Estimate + 1.96 * coef$`Std. Error`
    return(coef)
}


calc_ci_interval = function(model, x_grid){
    inv_func = family(model)$linkinv
    m_preds = predict(model, type = "link", se = T, newdata = x_grid) %>% 
        as_tibble() %>% 
        dplyr::mutate("upper_ci" = inv_func(fit + 1.96 * se.fit),
                      "lower_ci" = inv_func(fit - 1.96 * se.fit),
                      "fit" = inv_func(fit)) %>% 
        dplyr::select(fit, upper_ci, lower_ci)
    m_preds = cbind(x_grid, m_preds)
    return(m_preds)
}

calc_pred_interval <- function(model, x_grid) {
    data <- model$data
    n = nrow(data)
    set.seed(1)
    pred_values = purrr::map(1:10000, ~safe_calc_pred_interval_1boot(model, data, x_grid, n)) %>%
        purrr::transpose() %>% 
        purrr::pluck(1) %>% 
        do.call(rbind, .)
    boot_pi <- t(apply(pred_values, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
    colnames(boot_pi) = c("lower_pi", "upper_pi")
    pi_df = cbind(x_grid, boot_pi)
    return(pi_df)
}

calc_pred_interval_1boot = function(model, data, x_grid, n){
    boot = data[sample(seq_len(n), size = n, replace = TRUE), ]
    updated_model = update(model, data = boot)
    pred_lambda = predict(updated_model, type = "response", newdata = x_grid)
    pred_values = rpois(length(pred_lambda), lambda = pred_lambda)
    return(pred_values)
}
safe_calc_pred_interval_1boot = safely(calc_pred_interval_1boot)