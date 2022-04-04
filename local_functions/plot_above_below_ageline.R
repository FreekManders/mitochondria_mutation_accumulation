plot_above_below_ageline = function(higher_lower_muts){
    above_below_ageline_fig = ggplot(higher_lower_muts, aes(x = high_mut, y = nuclear_high_mut, fill = n)) +
        geom_raster() +
        geom_text(aes(label = round(n, 1)), size = 2) +
        labs(x = "Predicted mitochondrial mutation load", y = "Predicted nuclear mutation load", fill = "Number of samples") +
        theme_classic() +
        my_theme
    return(above_below_ageline_fig)
}
