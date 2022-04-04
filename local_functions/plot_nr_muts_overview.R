plot_nr_muts_overview = function(total_sample_mutation_freq){
    nr_muts_fig <- ggplot(total_sample_mutation_freq) +
        geom_histogram(aes(x = freq), 
                       fill = 'orange', 
                       binwidth = 1) +
        theme_BM(base_size = 30) +
        xlab('Mutations per sample') +
        ylab('Sample count') +
        facet_wrap(. ~ state)
    return(nr_muts_fig)
}

plot_nr_indels_overview = function(indel_df){
    indel_fig <- ggplot(indel_df) +
        geom_bar(aes(x = state, fill = patient), stat = 'count') +
        theme_BM(base_size = 30) +
        scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(length(unique(indel_df$patient))))
    return(indel_fig)
}
