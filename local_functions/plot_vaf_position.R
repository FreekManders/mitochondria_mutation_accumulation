plot_vaf_position = function(variant_df, muttype = "snv"){    
    
    if (muttype == "snv"){
        my_colors = c('C>A' = "#2EBAED", 'C>G' = "#000000", 'C>T' = "#DE1C14",'T>A' = '#D4D2D2',
                    'T>C' = "#ADCC54", 'T>G' = "#F0D0CE")

        # Get mutation type
        somatic_variant_gr = variant_df %>% 
            dplyr::mutate("chrom" = "MT") %>% 
            makeGRangesFromDataFrame(start.field = "mut_pos", end.field = "mut_pos", keep.extra.columns = T)
        variant_df$mut = factor(mut_type(somatic_variant_gr), levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
        
        legend_name = "Point mutation type"
    } else if (muttype == "indel"){
        variant_df$mut = ifelse(variant_df$indel == "insertion", "Insertion", "Deletion")
        legend_name = "Indel type"
        my_colors = c("Insertion" = "blue", "Deletion" = "red")
    }

    vaf_fig = ggplot(variant_df, 
                      aes(x = mut_pos, y = VAF)) +
        geom_point(size = 1, aes(col = mut)) +
        labs(x = 'mtDNA position (bp)',
             y = 'Variant Allele Frequency',
             colour = legend_name) +
        scale_color_manual(values=my_colors) +
        scale_x_continuous(limits=c(0,16550)) +
        scale_y_continuous(limits=c(0,1)) +
        theme_BM() +
        my_theme
    return(vaf_fig)
}
