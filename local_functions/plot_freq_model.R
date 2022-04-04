plot_freq_model = function(model, var){
    
    var = ensym(var)
    
    # Get data
    df = model$data
    
    # Calculate interval
    var_vals = eval(expr(`$`(df, !!var)))
    vars = expand.grid(var = unique(var_vals), "age" = seq(min(df$age), max(df$age), length.out = 20))
    colnames(vars)[1] = as.character(expr(!!var))
    ci_df = calc_ci_interval(glm_m, vars)
    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = !!var)) +
        geom_jitter(size = 2, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5), 
                   linetype = 2, 
                   alpha = 0.5) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1.5) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = !!var), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(fill = 'none') +
        annotate(geom='text', 
                 x = 0.19, 
                 y = 2,
                 label = paste('P',pval, sep = ' = '), 
                 size = 5) +
        theme_BM()
    return(ageline_fig)
}

plot_aml_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("state" = factor(levels(df$state), levels = levels(df$state)), "age" = seq(min(df$age), max(df$age), length.out = 20))
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)
    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = state)) +
        geom_jitter(size = 2, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5), 
                   linetype = 2, 
                   alpha = 0.5) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1.5) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = state), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(fill = 'none') +
        annotate(geom='text', 
                 x = 1, 
                 y = 3,
                 label = paste('P',pval, sep = ' = '), 
                 size = 5) +
        theme_BM()
    return(ageline_fig)
}

plot_pcawg_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("type" = factor(levels(df$type), levels = levels(df$type)), "age" = 0:90)
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 3)

    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = type)) +
        geom_jitter(size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5), 
                   linetype = 2, 
                   alpha = 0.25,
                   size = 0.25) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = type), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)", colour = "Type") +
        coord_cartesian(ylim = c(-0.5, 8)) +
        guides(fill = 'none') +
        scale_color_manual(values = c("Healthy blood" = "#E64B35FF",
                                      "Normal colon" = "#E64B35FF",
                                      "Normal liver" = "#E64B35FF",
                                      "Blood cancer" = "darkred",
                                      "Colon cancer" = "darkred",
                                      "Liver cancer" = "darkred"), limits = force) +
        scale_fill_manual(values = c("Healthy blood" = "#E64B35FF",
                                     "Normal colon" = "#E64B35FF",
                                     "Normal liver" = "#E64B35FF",
                                     "Blood cancer" = "darkred",
                                     "Colon cancer" = "darkred",
                                     "Liver cancer" = "darkred")) +
        annotate(geom='text', 
                 x = 12, 
                 y = 8,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        scale_x_continuous(breaks = c(0, 20, 40, 60, 80)) +
        scale_y_continuous(breaks = c(0:8)) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_hsct_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("state" = factor(levels(df$state), levels = levels(df$state)), "age" = 0:65)
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)

    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = state)) +
        geom_jitter(size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = state), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)", colour = "State") +
        guides(fill = 'none') +
        coord_cartesian(ylim = c(-0.5,6)) +
        scale_color_manual(values = c("Healthy blood" = "#E64B35FF",
                                      "HSCT recipient" = "#B09C85FF")) +
        scale_fill_manual(values = c("Healthy blood" = "#E64B35FF",
                                      "HSCT recipient" = "#B09C85FF")) +
        annotate(geom='text', 
                 x = 10, 
                 y = 6,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_cov_model2 = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("cov" = factor(levels(df$cov), levels = levels(df$cov)), "age" = 0:65)
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 3)

    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = cov)) +
        geom_jitter(size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = cov), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)", colour = "Sequencing depth") +
        guides(fill = 'none') +
        coord_cartesian(ylim = c(-0.5,6)) +
        annotate(geom='text', 
                 x = 10, 
                 y = 6,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_chemo_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("state" = factor(levels(df$state), levels = levels(df$state)), "age" = 0:65)
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)

    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = state)) +
        geom_jitter(size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = state), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)", colour = "State") +
        coord_cartesian(ylim = c(-0.5, 6)) +
        guides(fill = 'none') +
        scale_color_manual(values = c("Healthy blood" = "#E64B35FF", 
                                      "Diagnosis" = "#91D1C2FF",
                                      "HSPCs leukemia patients" = "#91D1C2FF",
                                      "Leukemia" = "#4DBBD5FF", 
                                      "Follow-up" = "#8491B4FF", 
                                      "Diagnosis 2" = "#00A087FF", 
                                      "2nd Leukemia" = "#3C5488FF"), limits = force) +
        scale_fill_manual(values = c("Healthy blood" = "#E64B35FF", 
                                      "Diagnosis" = "#91D1C2FF",
                                     "HSPCs leukemia patients" = "#91D1C2FF",
                                      "Leukemia" = "#4DBBD5FF", 
                                      "Follow-up" = "#8491B4FF", 
                                      "Diagnosis 2" = "#00A087FF", 
                                      "2nd Leukemia" = "#3C5488FF")) +
        annotate(geom='text', 
                 x = 10, 
                 y = 6,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    
    # Add extra pvalue if there are more than two state levels (healthy, chemo_DX and chemo_FU).
    if (length(levels(df$state)) >= 3){
        pval2 = calc_pval(model, 3)
        pval2 = round(as.numeric(pval2), 4)
        
        ageline_fig = ageline_fig +
            annotate(geom='text', 
                     x = 10, 
                     y = 5.5,
                     label = paste('P',pval2, sep = ' = '), 
                     size = 2)
    }
    return(ageline_fig)
}

plot_cb_chemo = function(cb_chemo_freq){
    aov_res = broom::tidy(aov(freq ~ chemo, data = cb_chemo_freq))
    pval = aov_res$p.value[1]
    
    cb_chemo_fig = ggplot(cb_chemo_freq, aes(x = chemo, y = freq)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(aes(colour = patient), size = 1, groupOnX = TRUE) +
        annotate("text", x = 5.5, y = 6, size = 2, label = paste0("P: ", round(pval, 4))) +
        scale_y_continuous(breaks = c(0, 2, 4, 6)) +
        theme_classic() +
        coord_cartesian(ylim = c(0,6)) +
        labs(x = "Treatment", y = "Base substitutions per clone") +
        guides(colour = 'none') +
        my_theme
    return(cb_chemo_fig)
}

plot_freq_cnv_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("cnv_mean" = seq(min(df$cnv_mean), max(df$cnv_mean), length.out = 20))
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)
    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = cnv_mean, y = freq)) +
        geom_jitter(aes(color = patient),
                    size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5), 
                   linetype = 2, 
                   alpha = 0.5) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "mtDNA Copy Number") +
        guides(color = 'none') +
        annotate(geom='text', 
                 x = 100, 
                 y = 3,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_recipient_time_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate fit
    df$fit = predict(model, level = 0)
    
    # Calculate pval
    pval = summary(model)$tTable[2, 5]
    pval <- formatC(as.numeric(pval),
                    digits = 4,
                    format = 'f')  
    
    # Create ageline figure
    fig <- ggplot(df, aes(x = time_transplant, y = cnv_mean)) +
        geom_jitter(aes(color = patient),
                    size = 1, 
                    width = 0.2, 
                    height = 0.1) +
        geom_line(aes(y = fit), 
                  size = 1) +
        labs(y = "mtDNA Copy Number", x = "Time after transplant (years)") +
        guides(fill = 'none', color = 'none') +
        annotate(geom='text', 
                 x = 0, 
                 y = 1750,
                 label = paste('P', pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    return(fig)
}


plot_bulk_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("bulk" = factor(levels(df$bulk), levels = levels(df$bulk)), "age" = seq(min(df$age), max(df$age), length.out = 20))
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)
    
    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = bulk)) +
        geom_jitter(size = 2, 
                    width = 0.2, 
                    height = 0.1) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1.5) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = bulk), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(fill = 'none') +
        annotate(geom='text', 
                 x = 1, 
                 y = 20,
                 label = paste('P',pval, sep = ' = '), 
                 size = 5) +
        theme_BM()
    return(ageline_fig)
}
calc_pval = function(model, i){
    pval <- formatC(summary(model)$coefficients[i,4],
                    digits = 4,
                    format = 'f')
    return(pval)
}

plot_fetus_model = function(model){
    
    df = model$data
    
    vars = expand.grid("state" = unique(df$state), "age" = seq(0.19, 0.268, by = 0.01))
    ci_df = calc_ci_interval(model, vars)

    pval <- calc_pval(model, 2)
    
    ageline_fig <- ggplot(df,
                           aes(x = age, y = freq, color = state)) +
        geom_jitter(size = 2, 
                    width = 0.001, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5), 
                   linetype = 2, 
                   alpha = 0.5) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1.5) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = state), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(fill = 'none') +
        annotate(geom='text', 
                 x = 0.2, 
                 y = 2,
                 label = paste('P',pval, sep = ' = '), 
                 size = 5) +
        theme_BM()
    return(ageline_fig)
}

plot_fetus_model2 = function(df){
 
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = state)) +
        geom_jitter(size = 1, 
                    width = 0.001, 
                    height = 0.1) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        labs(y = "Base substitutions per clone", x = "Age (years)", color = "State") +
        guides(fill = 'none') +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_ageline = function(model){
    df = model$data
    
    #ggeffects doesn't give difference between prediction and confidence for some reason.
    #Therefore values are calculated manually.
    x_grid = tibble("age" = 0:80)
    ci_df = calc_ci_interval(model, x_grid)
    pi_df = calc_pred_interval(model, x_grid)
    ci_pi_df = cbind(ci_df, pi_df[,-1])
    
    pval <- calc_pval(model, 2)
   
    ageline_fig <- ggplot(df,
                           aes(x = age, y = freq )) +
        geom_jitter(aes(col = patient), 
                    size = 1, 
                    width = 1, 
                    height = 0.25) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        geom_line(data = ci_pi_df, aes(y = fit), 
                  size = 1, 
                  color = 'red') +
        geom_ribbon(data = ci_pi_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit), alpha = 0.2) +
        geom_ribbon(data = ci_pi_df, aes(ymin = lower_pi, ymax = upper_pi, y = fit), alpha = 0.1) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(color = 'none') +
        coord_cartesian(ylim = c(-0.5, 7)) +
        annotate(geom='text', 
                 x = 12, 
                 y = 7,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_liver_ageline = function(df){
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq )) +
        geom_jitter(aes(col = patient), 
                    size = 1, 
                    width = 1, 
                    height = 0.25) +
        geom_hline(yintercept = c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5), 
                   linetype = 2, 
                   alpha = 0.5,
                   size = 0.25) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(color = 'none') +
        coord_cartesian(ylim = c(-0.5, 7), xlim = c(0, 65)) +
        scale_color_manual(values = c("cornflowerblue", "goldenrod2")) +
        theme_BM() +
        my_theme
    return(ageline_fig)
}

plot_depth_model = function(model){
    
    # Get data
    df = model$data
    
    # Calculate interval
    vars = expand.grid("cov" = factor(levels(df$cov), levels = levels(df$cov)), "age" = seq(0.19, 0.268, by = 0.01))
    ci_df = calc_ci_interval(model, vars)
    
    # Calculate pval
    pval = calc_pval(model, 2)
    pval <- formatC(as.numeric(pval),
                    digits = 4,
                    format = 'f')  
    
    # Create ageline figure
    ageline_fig <- ggplot(df,
                          aes(x = age, y = freq, color = cov)) +
        geom_jitter(size = 2, 
                    width = 0.001, 
                    height = 0.1) +
        geom_line(data = ci_df, aes(y = fit), 
                  size = 1.5) +
        geom_ribbon(data = ci_df, aes(ymin = lower_ci, ymax = upper_ci, y = fit, fill = cov), color = NA, alpha = 0.2) +
        labs(y = "Base substitutions per clone", x = "Age (years)") +
        guides(fill = 'none') +
        annotate(geom='text', 
                 x = 0.2, 
                 y = 2,
                 label = paste('P',pval, sep = ' = '), 
                 size = 5) +
        theme_BM()
    return(ageline_fig)
}
