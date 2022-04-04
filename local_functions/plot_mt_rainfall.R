plot_mt_rainfall = function(variants){
    COLORS6 = c('C>A' = "#2EBAED", 'C>G' = "#000000", 'C>T' = "#DE1C14",'T>A' = '#D4D2D2',
                'T>C' = "#ADCC54", 'T>G' = "#F0D0CE")
    
    gr = to_granges_mt(variants)
    healthy.colors <- getVariantsColors(gr$ref, gr$alt, color.table = COLORS6)
    
    kp <- plotKaryotype(chromosomes = 'chrM', 
                        genome = 'hg38', 
                        ideogram.plotter = NULL,
                        labels.plotter = NULL)
    kpAddCytobandsAsLine(kp)
    kpAddChromosomeNames(kp, srt=45)
    kpAxis(kp, ymax = 7, 
           tick.pos = 1:7, 
           cex = 2)
    kpPlotRainfall(kp, 
                   data = gr, 
                   col = healthy.colors, 
                   cex = 1.2)
    kpAddLabels(kp, labels = c("Distance between mutations (log10)"), 
                srt=90, 
                pos=1, 
                label.margin = 0.06, 
                cex = 1.5)
    invisible(0)
}

plot_save_rainfall = function(variants, name){
    pdf(file.path(plotdir, "rainfall", paste0("rainfall_", name, ".pdf")))
    plot_mt_rainfall(variants)
    dev.off()
}