# Plot average MT coverage
plot_cov <- function(df, include_pat, dir) {
    plotname = paste('mitochondrial_coverage', include_pat, sep = '_')
    df_patient = dplyr::filter(df, patient %in% include_pat)
    df_patient_l = split(df_patient, df_patient$sample)
    
    # Calculate mean coverage of windows
    df_patient = purrr::map(df_patient_l, calc_coverage_windows, window_size = 20) %>% 
        bind_rows()
    
    max_y = max(df_patient$coverage)
    
    plot = ggplot(df_patient, aes(x = position, y = coverage)) +
        geom_point(aes(col = sample), pch = 20, size = 1, alpha = 0.2) +
        labs(x = 'mtDNA position (bp)',
             y= 'Mean MT Coverage',
             colour = "Sample") +
        coord_cartesian(ylim = c(0, max_y)) +
        theme_BM() +
        my_theme
    ggsave(paste(plotname, 'pdf', sep = '.'), plot, device = 'pdf', path = dir)
    saveRDS(plot, file.path(dir, paste0(plotname, ".rds")))
    return(plot)
}

