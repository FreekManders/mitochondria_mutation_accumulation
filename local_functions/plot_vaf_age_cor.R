plot_vaf_age_cor = function(variants){
    
    # Get variants with max vaf in that patient.
    vaf_tb = variants %>% 
        dplyr::group_by(patient) %>% 
        dplyr::mutate(max_vaf = VAF == max(VAF))
    max_vaf_tb = dplyr::filter(vaf_tb, max_vaf & bulk == "clone")
    
    r2 = (cor(max_vaf_tb$VAF, max_vaf_tb$age))^2
    r2 <- formatC(r2,
                    digits = 4,
                    format = 'f')
    max_vaf_m_tb = broom::tidy(lm(VAF ~ age, data = max_vaf_tb))
    pval <- formatC(max_vaf_m_tb$p.value[2],
                    digits = 4,
                    format = 'f')
    
    maxvaf_cor_fig = plot_single_vaf_age_cor(max_vaf_tb, r2, pval)
    return(maxvaf_cor_fig)
}
plot_single_vaf_age_cor = function(vaf_tb, r2, pval){
    fig = ggplot(vaf_tb, aes(x = age, y = VAF)) +
        geom_jitter(width = 1, size = 1) +
        geom_smooth(method = "lm", color = "black", size = 1) +
        annotate("text", x = 2, y = 0.9, label = paste0("R2: ", r2), size = 2) +
        annotate("text", x = 2, y = 0.95, label = paste0("P: ", pval), size = 2) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(y = "Maximum VAF", x = "Age (years)") +
        theme_BM() +
        my_theme
    return(fig)
}
