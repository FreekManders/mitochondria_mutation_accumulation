plot_mean_hsct_freq = function(hsct_freq){
    
    # Remove non hsct and calculate mean
    mean_hsct_freq = hsct_freq %>% 
        dplyr::filter(hsct != "No") %>% 
        dplyr::group_by(patient, hsct) %>% 
        dplyr::summarise(mean_freq = mean(freq), .groups = "drop")
    
    # plot figure
    mean_hsct_fig = ggplot(mean_hsct_freq, aes(x = patient, y = mean_freq, fill = hsct)) +
        geom_bar(position = "dodge", stat = "identity") +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(mean_hsct_fig)
}


plot_mean_dx2_freq = function(freq){
    
    # Remove non chemo samples
    freq = freq %>% 
        dplyr::filter(chemo == "DX1" | chemo == "DX2")
    
    # plot figure
    mean_chemo_fig = ggplot(freq, aes(x = patient, colour = chemo, y = freq)) + 
        geom_boxplot(outlier.shape = NA) + 
        geom_quasirandom(groupOnX = T, dodge.width = 0.75) +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(mean_chemo_fig)
}


plot_mean_chemo_freq = function(freq){
    
    # Remove non chemo samples
    freq = freq %>% 
        dplyr::filter(chemo == "DX1" | chemo == "FU")
    
    # plot figure
    mean_chemo_fig = ggplot(freq, aes(x = patient, colour = chemo, y = freq)) + 
        geom_boxplot(outlier.shape = NA) + 
        geom_quasirandom(groupOnX = T, dodge.width = 0.75) +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(mean_chemo_fig)
}

plot_mean_chemo_cnv = function(cnv){
    
    # Remove non chemo samples
    cnv = cnv %>% 
        dplyr::filter(chemo == "chemo_DX" | chemo == "chemo_FU")
    
    # plot figure
    mean_chemo_fig = ggplot(cnv, aes(x = patient, colour = chemo, y = cnv_mean)) + 
        geom_boxplot(outlier.shape = NA, lwd = 0.25) + 
        geom_quasirandom(groupOnX = T, dodge.width = 0.75) +
        theme_classic() +
        theme(text = element_text(size = 20))
    return(mean_chemo_fig)
}

plot_mean_hsct_cnv = function(hsct_cnv){
    
    # Remove non hsct and calculate mean
    hsct_cnv = hsct_cnv %>% 
        dplyr::filter(patient %in% c("HSCT2", "HSCT3", "HSCT4", "HSCT13", "HSCT14")) %>%
        dplyr::mutate(patient = factor(patient, levels = c("HSCT2", "HSCT3", "HSCT4", "HSCT13", "HSCT14")),
                      hsct = dplyr::recode(hsct, "donor" = "Donor", "recipient" = "Recipient")) %>% 
        dplyr::group_by(patient, hsct)
    

    # plot figure
    fig = ggplot(hsct_cnv, aes(x = hsct, y = cnv_mean)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_point() +
        facet_grid(. ~ patient) +
        labs(y = "mtDNA Copy Number", x = "") +
        theme_BM() +
        my_theme +
        theme(strip.background = element_rect(size = 0.5, fill = NA))
    return(fig)
}
