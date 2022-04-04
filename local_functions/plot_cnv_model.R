plot_cnv_ageline = function(model){
    
    # Get data
    df = model$data
    
    
    
    # Calculate fit and confidence interval
    ci_df = ggpredict(model, "age [0:80]")
    colnames(ci_df)[1] = "age"
    
    # Use in case lme model gives issues.
    m2 = lmer(cnv_mean ~ age + (0 + age | patient),
        data = model$data)
    
    pi_df = ggpredict(m2, "age [0:80]", type = "re")
    colnames(pi_df)[1] = "age"
    
    #Get pvalue
    sum_m <- summary(model)
    pval <- formatC(sum_m$tTable[2,5],
                    digits = 4,
                    format = 'f')
    
    # Plot figure
    fig_healthy_cnv = ggplot(df,
                             aes(x = age, 
                                 y = cnv_mean )) +
        geom_jitter(aes(col = patient),
                    size = 1, 
                    width = 1, 
                    height = 0.3) +
        geom_line(data = ci_df,
                  aes(y = predicted), 
                  color = 'red', 
                  size = 1) +
        geom_ribbon(data = ci_df, aes(ymin = conf.low, ymax = conf.high, y = predicted), color = NA, alpha = 0.2) +
        geom_ribbon(data = pi_df, aes(ymin = conf.low, ymax = conf.high, y = predicted), color = NA, alpha = 0.1) +
        annotate(geom='text', 
                 x = 10, 
                 y = 1600,
                 label = paste('P',pval, sep = ' = '), 
                 size = 2) +
        coord_cartesian(ylim = c(0, 1600), xlim = c(0, 80)) +
        ylab("mtDNA Copy Number") +
        xlab("Age (years)") +
        theme_BM() +
        guides(color = 'none') +
        my_theme
    return(fig_healthy_cnv)
}


plot_cnv_liver_ageline = function(df){
    fig_healthy_cnv = ggplot(df,
                             aes(x = age, 
                                 y = cnv_mean )) +
        geom_jitter(aes(col = patient),
                    size = 1, 
                    width = 1, 
                    height = 0.3) +
        coord_cartesian(ylim = c(0, 1600), xlim = c(0, 65)) +
        ylab("mtDNA Copy Number") +
        xlab("Age (years)") +
        scale_color_manual(values = c("cornflowerblue", "goldenrod2")) +
        theme_BM() +
        guides(color = 'none') +
        my_theme
    return(fig_healthy_cnv)
}

plot_cnv_model = function(model, var = "state", col = "patient", plot_p = TRUE, y_annotate = 1500, remove_guide = TRUE, size = 1){
    
    # Make symbol of var
    var = ensym(var)
    col = ensym(col)
    
    # Get data
    df = model$data
    
    #Get pvalue
    sum_m <- summary(model)
    pval <- formatC(sum_m$tTable[2,5],
                    digits = 4,
                    format = 'f')
    
    # Plot figure
    fig = ggplot(df, aes(x = !!var, y = cnv_mean)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(aes(col = !!col),
                         size = size, 
                         width = 0.3) +
        ylab("mtDNA Copy Number") +
        xlab("") +
        theme_BM() +
        my_theme
    
    if (remove_guide){
        fig = fig +
            guides(color = 'none')
    }
    
    
    if (plot_p == TRUE){
        fig = fig +
            annotate(geom='text', 
                     x = 1.5, 
                     y = y_annotate,
                     label = paste('P',pval, sep = ' = '), 
                     size = 2)
    }
    return(fig)
}

plot_aml_cnv_model = function(df){
    
    #Get pvalue
    res = broom::tidy(wilcox.test(cnv_mean ~ source, data = aml_cnv))
    pval <- formatC(res$p.value,
                    digits = 4,
                    format = 'f')
    
    # Plot figure
    fig = ggplot(df, aes(x = source, y = cnv_mean)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(aes(col = patient),
                         size = 1, 
                         width = 0.3) +
        ylab("mtDNA Copy Number") +
        xlab("Sample source") +
        theme_BM() +
        guides(color = 'none') +
        my_theme +
        annotate(geom='text', 
                     x = 1.5, 
                     y = 2700,
                     label = paste('P',pval, sep = ' = '), 
                     size = 2)
    
    return(fig)
}


plot_cb_chemo_cnv = function(cb_chemo_cnv){
    aov_res = broom::tidy(aov(cnv_mean ~ chemo, data = cb_chemo_cnv))
    pval = aov_res$p.value[1]

    cb_chemo_fig = ggplot(cb_chemo_cnv, aes(x = chemo, y = cnv_mean)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(aes(colour = patient), size = 1, groupOnX = TRUE) +
        annotate("text", x = 5.5, y = 1350, size = 2, label = paste0("P: ", round(pval, 4))) +
        labs(x = "Treatment", y = "mtDNA Copy Number") +
        guides(colour = 'none') +
        theme_classic() +
        my_theme
    return(cb_chemo_fig)
}
