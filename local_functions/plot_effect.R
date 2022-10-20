plot_effect = function(annotate_df){

    # Set label order
    annotate_df = dplyr::filter(annotate_df, bulk == "clone" & !is.na(effect))
    annotate_df$effect <- factor(annotate_df$effect, levels = c('LOW','MODERATE','HIGH'))
    
    plot <- ggplot(annotate_df, aes(x = effect)) +
        geom_bar(stat = 'count', aes(fill = effect)) +
        labs(x = 'Predicted effect', y = "Nr. substitutions") +
        scale_fill_manual(values = c('green','orange','red'), drop = FALSE, guide = "none") +
        scale_x_discrete(drop = FALSE) +
        facet_wrap(. ~ state_name, nrow =1) +
        theme_BM() +
        my_theme +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5, size = 6, family = "Arial"),
              panel.spacing = unit(0.1, "lines"),
              strip.text.x = element_text(size = 3, family = "Arial"))
    return(plot)
}

plot_gene = function(annotate_df, length_correction = FALSE){
    
    # Change NA to None
    annotate_df$gene[is.na(annotate_df$gene)] <- c('None')
    
    # Count mutations per gene
    count_tb = annotate_df %>%
        dplyr::filter(bulk == "clone") %>%
        dplyr::filter(!duplicated(variant)) %>% # Remove duplicate variants. Because these are shared within a single donor, the duplicates can be removed without causing issues with the groups..
        dplyr::mutate(strand = dplyr::recode(strand, "1" = "+", "-1" = "-"),
                      gene = paste0(gene, " (", strand, ")"),
                      gene = factor(gene, levels = sort(unique(gene))),
                      state_name = droplevels(state_name)) %>% 
        dplyr::group_by(state_name, gene, length, .drop = FALSE) %>% 
        dplyr::count(.drop = FALSE) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(gene = ifelse(gene == "None (NA)", "Not genic", as.character(gene)))
    
    
    # Set plot label
    fill_label = "Nr. mutations"
    if(length_correction){
        count_tb = dplyr::mutate(count_tb, n = ifelse(is.na(length), n, 1000 * n / length))
        #CHANGE LABEL n
        fill_label = "Nr. substitutions / gene length (kb)"
    }
    
    # Create plot
    fig = ggplot(count_tb, aes(x = state_name, y = gene, fill = n)) +
        geom_raster() +
        geom_text(aes(label = round(n, 1)), size = 2) +
        scale_fill_distiller(palette = "RdYlBu") +
        labs(fill = fill_label, x = "", y = "Gene") +
        theme_BM() +
        my_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6, family = "Arial"))
    
    return(fig)
}
