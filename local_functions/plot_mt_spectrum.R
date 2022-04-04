plot_mt_indel_spectrum = function(variants, ref_genome){
    gr = to_granges_mt(variants)
    gr = get_indel_context(gr, ref_genome)
    counts = count_indel_contexts(gr)
    fig = plot_indel_contexts(counts, condensed = TRUE)
    return(fig)
}

plot_mt_96_profile = function(variants, ref_genome, name){
    # Subset to C/T and G/A
    variants_c_t <- subset(variants, ref %in% c('C','T'))
    variants_g_a <- subset(variants, ref %in% c('G','A'))
    
    # Transform into GRangesList
    gr = to_granges_mt(variants)
    ct_gr = to_granges_mt(variants_c_t)
    ga_gr = to_granges_mt(variants_g_a)
    
    grl = GRangesList(gr, ct_gr, ga_gr)
    names(grl) = paste0(c("", "ct_", "ga_"), name)
    
    # Create matrix
    mut_mat = mut_matrix(grl, ref_genome)
    
    # Create figure
    fig = plot_96_profile(mut_mat, condensed = TRUE)
    return(fig)
}

plot_save_spectra = function(variants, name){
    spectrum_fig = plot_mt_spectrum(variants, ref_genome, name = name)
    ggsave(file.path(plotdir, "spectra", paste0("spectrum_", name, ".pdf")), spectrum_fig)
    
    profile_fig = plot_mt_96_profile(variants, ref_genome, name = name)
    ggsave(file.path(plotdir, "spectra", paste0("profile96_", name, ".pdf")), profile_fig)
}

plot_mt_spectrum = function(variants, ref_genome, name, remove_duplis = TRUE){
    
    # Get split spectra
    spec_l = get_mt_split_spectra(variants, ref_genome, remove_duplis)
    
    # plot combined spectrum
    spec_combined <- rbind(spec_l$spec, spec_l$spec_c_t, spec_l$spec_g_a)
    exp_names = paste0(name, c("", ": Light strand", ": Heavy strand"))
    spectrum_fig <- plot_spectrum(spec_combined, CT = TRUE, by = exp_names, error_bars = "none")
    return(spectrum_fig)
}

get_mt_split_spectra = function(variants, ref_genome, remove_duplis){
    # Subset to C/T and G/A
    variants_c_t <- subset(variants, ref %in% c('C','T'))
    variants_g_a <- subset(variants, ref %in% c('G','A'))
    
    # Transform into GRanges
    gr = to_granges_mt(variants, remove_duplis)
    ct_gr = to_granges_mt(variants_c_t, remove_duplis)
    ga_gr = to_granges_mt(variants_g_a, remove_duplis)
    
    # Get single type occurrence
    spec <- single_type_occurence(gr, ref_genome)
    spec_c_t <- single_type_occurence(ct_gr, ref_genome)
    spec_g_a <- single_type_occurence(ga_gr, ref_genome)
    
    spec_l = list("spec" = spec, "spec_c_t" = spec_c_t, "spec_g_a" = spec_g_a)
    return(spec_l)
}

to_granges_mt = function(variants, remove_duplis = TRUE){
    
    if (remove_duplis){
        # Remove duplicate variants. Unique can't be used later, because
        # this doesn't take ref and alt into account.
        dupli_f = duplicated(variants[,c("mut_pos", "ref", "alt")])
        variants = variants[!dupli_f,]
    }
    gr <- GRanges(seqnames = rep("chrM", nrow(variants)),
                  ranges = IRanges(start = variants$mut_pos, end = variants$mut_pos),
                  strand = rep('*', nrow(variants)),
                  ref = variants$ref,
                  alt = variants$alt)
    GenomeInfoDb::genome(gr) = "hg38"
    return(gr)
}
