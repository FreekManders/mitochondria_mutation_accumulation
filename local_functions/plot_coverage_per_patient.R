plot_coverage_per_patient = function(cnv_df){
    
    patient_states = cnv_df %>% dplyr::filter(bulk == "clone") %>% dplyr::group_by(patient) %>% dplyr::summarise(state_name = state_name[[1]])
    patient_order = patient_states[order(patient_states$state_name),"patient", drop = TRUE]
    
    #cnv_df$patient = factor(cnv_df$patient, levels = unique(cnv_df$patient))
    cnv_df$patient = factor(cnv_df$patient, levels = patient_order)
    
    coverage_fig = ggplot(cnv_df, 
                          aes(x = patient, 
                              y = mt_mean)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(aes(col = state_name), 
                         size = 0.5) +
        geom_hline(yintercept = 1000, 
                   col = 'red', 
                   linetype = 2) +
        labs(x = "Donor", y = "Mean MT Coverage (log10)", col = "State") +
        theme_BM() +
        scale_y_log10(labels = scales::comma) +
        my_theme +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.key.size = unit(0.3, 'cm'))
    return(coverage_fig)
}
