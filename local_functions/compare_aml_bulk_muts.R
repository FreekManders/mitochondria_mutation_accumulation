compare_aml_bulk_muts = function(target_muts){
    
    # Calculate difference between clone and bulk
    diff_tb = target_muts %>% 
        dplyr::select(patient, bulk, freq) %>% 
        tidyr::pivot_wider(names_from = "bulk", values_from = "freq") %>% 
        dplyr::mutate(diff = clone - bulk)
    
    # Do statistical test
    wilcox_res = broom::tidy(wilcox.test(diff_tb$clone, diff_tb$bulk, paired = TRUE, exact = FALSE))
    
    # Create figure
    fig = ggplot(diff_tb, aes(y = diff, x = "")) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom() +
        annotate("text", y = 2, x = 0.5, label = paste0("P: " ,round(wilcox_res$p.value, 3))) +
        labs(y = "AML muts - bulk") +
        theme_classic()
    return(fig)
}


compare_aml_bulk_cnv = function(target_cnv){
    
    # Calculate difference between clone and bulk
    diff_tb = target_cnv %>% 
        dplyr::select(patient, bulk, cnv_mean) %>% 
        tidyr::pivot_wider(names_from = "bulk", values_from = "cnv_mean") %>% 
        dplyr::mutate(diff = clone - bulk)
    
    # Do statistical test
    wilcox_res = broom::tidy(wilcox.test(diff_tb$clone, diff_tb$bulk, paired = TRUE, exact = FALSE))
    
    target_cnv_2plot = target_cnv %>% 
        dplyr::mutate(bulk = ifelse(bulk == "bulk", "Bulk blood", "AML"),
                      bulk = factor(bulk, levels = c("Bulk blood", "AML")))
    # Create figure
    fig = ggplot(target_cnv_2plot, aes(x = bulk, y = cnv_mean)) +
        geom_point(size = 1) +
        geom_line(aes(group = patient)) +
        annotate("text", y = 500, x = 1.5, label = paste0("P: " ,round(wilcox_res$p.value, 4)), size = 2) +
        labs(y = "mtDNA Copy Number", x = "") +
        theme_classic() +
        my_theme

    return(fig)
}

compare_aml_bulk_mt_reads = function(target_cnv){
    
    # Calculate difference between clone and bulk
    diff_tb = target_cnv %>% 
        dplyr::select(patient, bulk, mt_mean) %>% 
        tidyr::pivot_wider(names_from = "bulk", values_from = "mt_mean") %>% 
        dplyr::mutate(diff = clone - bulk)
    
    # Do statistical test
    wilcox_res = broom::tidy(wilcox.test(diff_tb$clone, diff_tb$bulk, paired = TRUE, exact = FALSE))
    
    target_cnv_2plot = target_cnv %>% 
        dplyr::mutate(bulk = ifelse(bulk == "bulk", "Bulk blood", "AML"),
                      bulk = factor(bulk, levels = c("Bulk blood", "AML")))
    # Create figure
    fig = ggplot(target_cnv_2plot, aes(x = bulk, y = mt_mean)) +
        geom_point(size = 1) +
        geom_line(aes(group = patient)) +
        annotate("text", y = 12000, x = 1.5, label = paste0("P: " ,round(wilcox_res$p.value, 5)), size = 2) +
        labs(y = "Mean MT Coverage", x = "") +
        theme_classic() +
        my_theme
    
    return(fig)
}
