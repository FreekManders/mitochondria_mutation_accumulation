mito_cossim_spectra = function(variants1, variants2, remove_duplis1, remove_duplis2){
    
    # Get spectra
    spec1_l = get_mt_split_spectra(variants1, ref_genome, remove_duplis1)
    spec2_l = get_mt_split_spectra(variants2, ref_genome, remove_duplis2)
    
    # Calculate cosine similarities
    cossim = cos_sim(unlist(spec1_l$spec), unlist(spec2_l$spec))
    cossim_c_t = cos_sim(unlist(spec1_l$spec_c_t), unlist(spec2_l$spec_c_t))
    cossim_g_a = cos_sim(unlist(spec1_l$spec_g_a), unlist(spec2_l$spec_g_a))
    cossim_l = list("cossim" = cossim,
                    "cossim_c_t" = cossim_c_t,
                    "cossim_g_a" = cossim_g_a)
    return(cossim_l)
}

plot_save_pcawg_spectra_tissue = function(pcawg_snv, variants, name, cancer_name, high_mut_samples, cancer_types){
    plot_save_pcawg_spectra(pcawg_snv, variants, "pcawg", name, type = "all_pcawg")
    
    blood_pcawg_snv = pcawg_snv %>% 
        dplyr::filter(cancer_type %in% cancer_types)
    plot_save_pcawg_spectra(blood_pcawg_snv, variants, cancer_name, name, type = "tissue_pcawg")
    
    high_mut_snv = blood_pcawg_snv %>% 
        dplyr::filter(sample_id %in% high_mut_samples)
    plot_save_pcawg_spectra(high_mut_snv, variants, paste0(cancer_name, "_high_mut"), name, type = "high_mut_pcwag")
}

plot_save_pcawg_spectra = function(pcawg_variants, healthy_variants, name1, name2, type){
    
    out_dir = file.path(plotdir, "pcawg", name2, paste0(type, "_spectra"))
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # pcawg spectrum
    pcawg_spectrum_fig = plot_mt_spectrum(pcawg_variants, ref_genome, name1, remove_duplis = FALSE) +
        spectrum_theme
    ggsave(paste0(name1, "_spectrum.pdf"), pcawg_spectrum_fig)
    saveRDS(pcawg_spectrum_fig, paste0(name1, "_spectrum.rds"))
    
    # pcawg compare spectra
    mito_cossims = mito_cossim_spectra(healthy_variants, pcawg_variants, 
                                       remove_duplis1 = TRUE, remove_duplis2 = FALSE)
    healthy_spectrum_fig = plot_mt_spectrum(healthy_variants, ref_genome, name2, remove_duplis = TRUE) +
        spectrum_theme
    saveRDS(healthy_spectrum_fig, "healthy_spectrum.rds")
    pcawg_spectrum_comparison_fig = ggarrange(healthy_spectrum_fig, 
                                              pcawg_spectrum_fig, 
                                              nrow = 2, 
                                              common.legend = TRUE)
    ggsave(paste0(name1, "_spectrum_comparison.pdf"), pcawg_spectrum_comparison_fig)
    saveRDS(pcawg_spectrum_comparison_fig, paste0(name1, "_spectrum_comparison.rds"))
    anno_pcawg_spectrum_comparison_fig = pcawg_spectrum_comparison_fig +
        annotate("text", x = 0.5, y = 0.472, size = 6, label = paste0("cossim: ", 
                                                                      round(mito_cossims$cossim, 4), 
                                                                      "; ", 
                                                                      round(mito_cossims$cossim_c_t, 4),
                                                                      "; ",
                                                                      round(mito_cossims$cossim_g_a, 4)))
    ggsave(paste0(name1, "_spectrum_comparison_anno.pdf"), anno_pcawg_spectrum_comparison_fig)
    
    
    # pcawg profile comparison
    healthy_gr = to_granges_mt(healthy_variants, remove_duplis = TRUE)
    pcawg_gr = to_granges_mt(pcawg_variants, remove_duplis = FALSE)
    
    grl = GRangesList(name2 = healthy_gr, name1 = pcawg_gr)
    mut_mat = mut_matrix(grl, ref_genome)
    pcawg_profile_comparison_fig = plot_compare_profiles(mut_mat[,1], mut_mat[,2], profile_names = names(grl))
    ggsave(paste0(name1, "_profile_comparison.pdf"), pcawg_profile_comparison_fig)
    
    invisible(0)
}



compare_pcawg_freq = function(pcawg_freq, state, name, total_sample_mutation_freq, pcawg_cancer_types, cancer_name, healthy_glm, pcawg){
    
    # Compare healthy to tissue pcawg data including the age variable.
    healthy_age_freq = total_sample_mutation_freq %>% 
        dplyr::filter(state == !!state & !(patient == 'AHH1') & bulk == 'clone') %>% 
        dplyr::mutate(type = !!name) %>% 
        dplyr::select(type, freq, age, sample)
    
    
    if (state == "healthy"){
        healthy_age_freq = total_sample_mutation_freq %>% 
            dplyr::filter((state == "healthy" | chemo %in% c("DX1_cancer")) & !(patient == 'AHH1') & bulk == 'clone') %>% 
            dplyr::mutate(type = dplyr::recode(state, "AML" = "Blood cancer", "ALL" = "Blood cancer", "healthy" = "Healthy blood")) %>% 
            dplyr::select(type, freq, age, sample)
    }
    
    pcawg_age_freq = pcawg %>%
        dplyr::select(type = histology_abbreviation, freq = nr_snvs, age = donor_age_at_diagnosis, sample = sample_id) %>% 
        dplyr::mutate(type = ifelse(type %in% pcawg_cancer_types, 
                                    cancer_name,
                                    as.character(type))) %>% 
        dplyr::filter(type %in% c(name, cancer_name) & age != "UNK") %>% 
        dplyr::mutate(age = as.numeric(age),
                      age = age + 0.75) #Count from conception
    
    # Combine healthy with pcawg data
    cancer_healthy_age_freq = rbind(pcawg_age_freq, healthy_age_freq) %>% 
        dplyr::mutate(type = fct_relevel(type, name))
    
    glm_m = tryCatch(expr = {glm(freq ~ age + type, data = cancer_healthy_age_freq, family = poisson(link = "identity"))},
             error = function(e){
                 glm(freq ~ age + type, data = cancer_healthy_age_freq, family = poisson(link = "identity"), start = c(coefficients(healthy_glm), 0))
             })
             
    
    pcawg_model_fig = plot_pcawg_model(glm_m)
    ggsave(file.path(plotdir, "pcawg", name, paste0(cancer_name, "pcawg_age_model.pdf")), pcawg_model_fig)
    saveRDS(pcawg_model_fig, file.path(plotdir, "pcawg", name, paste0(cancer_name, "pcawg_age_model_fig.rds")))
    saveRDS(glm_m, file.path(plotdir, "pcawg", name, paste0(cancer_name, "pcawg_age_model.rds")))
    
    # Check if high mut samples have a different histology
    high_mut_samples = pcawg_age_freq$sample[pcawg_age_freq$freq >= 5]
    histo_high_low = pcawg %>% 
        dplyr::filter(histology_abbreviation %in% pcawg_cancer_types) %>% 
        dplyr::mutate(high_mut = ifelse(sample_id %in% high_mut_samples, "High mut", "Low mut"),
                      high_mut = factor(high_mut, levels = unique(high_mut)),
                      histology_abbreviation = factor(histology_abbreviation, levels = unique(histology_abbreviation))) %>%
        dplyr::select(histology_abbreviation, high_mut)
    histo_high_low_m = histo_high_low %>%
        dplyr::count(histology_abbreviation, high_mut, .drop = FALSE) %>% 
        tidyr::pivot_wider(names_from = "histology_abbreviation", values_from = "n") %>% 
        as.data.frame() %>%
        tibble::column_to_rownames("high_mut") %>% 
        as.matrix()
    chisq_res = histo_high_low_m %>% 
        chisq.test(simulate.p.value = TRUE) %>% 
        broom::tidy()
    write_tsv(chisq_res, file.path(r_wd, paste0(name, "_high_mut_samples_histology_test.txt")))
    
    
    return(high_mut_samples)
}

compare_pcawg_cnv = function(pcawg_copy_number, cnv_df, state, name, cancer_name, pcawg_cancer_types){
    enquo(state)
    enquo(name)
    
    healthy_cnv = cnv_df %>% 
        dplyr::filter(bulk == "clone" & state == !!state) %>%
        dplyr::mutate(type = !!name) %>% 
        dplyr::select(type, copy_number = cnv_mean, patient)
    
    
    
    if (state == "healthy"){
        healthy_cnv = cnv_df %>% 
            dplyr::filter(bulk == "clone" & (state == "healthy" | chemo %in% c("DX1_cancer"))) %>%
            dplyr::mutate(type = dplyr::recode(state, "AML" = "Blood cancer", "ALL" = "Blood cancer", "healthy" = "Healthy blood")) %>% 
            dplyr::select(type, copy_number = cnv_mean, patient)
    }
    
    pcawg_healthy_cnv = rbind(pcawg_copy_number, healthy_cnv) %>% 
        dplyr::mutate(type = fct_relevel(type, name))
    
    cnv_per_cancer_fig = ggplot(pcawg_healthy_cnv, aes(x = type, y = copy_number, fill = type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(groupOnX = TRUE, size = 0.5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90))
    ggsave(file.path(plotdir, "pcawg", name, "pcawg_cnv_per_cancer.pdf"), cnv_per_cancer_fig)
    
    # Cancers are significantly different from each other. (Not new.)
    summary(aov(copy_number ~ type, data = pcawg_healthy_cnv))
    
    # Compare healthy to all pcawg cnv
    pcawg_vs_healthy_cnv = dplyr::mutate(pcawg_healthy_cnv, 
                                         type = ifelse(type == !!name, !!name, "pcawg"))
    
    pval_df = broom::tidy(wilcox.test(copy_number ~ type, data = pcawg_vs_healthy_cnv))
    pcawg_vs_healthy_cnv_fig = ggplot(pcawg_vs_healthy_cnv, aes(x = type, y = copy_number, fill = type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(groupOnX = TRUE, size = 0.5) +
        annotate("text", x = 1, y = 1500, label = paste0("P: ", round(pval_df$p.value, 4))) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90))
    ggsave(file.path(plotdir, "pcawg", name, "pcawg_vs_healthy_cnv.pdf"), pcawg_vs_healthy_cnv_fig)
    
    
    
    # Compare healthy cnv to tissue pcawg. Also compare their variance
    tissue_pcawg_vs_healthy_cnv = pcawg_healthy_cnv %>% 
        dplyr::mutate(type = ifelse(type %in% c(pcawg_cancer_types), 
                                    !!cancer_name, 
                                    as.character(type))) %>% 
        dplyr::filter(type %in% c(name, cancer_name)) %>% 
        dplyr::mutate(type = fct_relevel(type, name))
    
    m = lme(copy_number ~ type, random =  ~1 | patient, data = tissue_pcawg_vs_healthy_cnv)
    pval = summary(m)$tTable[2, "p-value"]
    pval <- formatC(pval,
                    digits = 4,
                    format = 'f')
    
    #pval_df = broom::tidy(wilcox.test(copy_number ~ type, data = blood_pcawg_vs_healthy_cnv))
    pval_df2 = broom::tidy(fligner.test(copy_number ~ type, data = tissue_pcawg_vs_healthy_cnv))
    tissue_pcawg_vs_healthy_cnv_fig = ggplot(tissue_pcawg_vs_healthy_cnv, aes(x = type, y = copy_number, fill = type)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(groupOnX = TRUE, size = 0.5) +
        annotate("text", x = 1.5, y = 1700, label = paste0("P: ", pval), size = 2) +
        labs(fill = "Type", x = "", y = "mtDNA Copy Number") +
        coord_cartesian(ylim = c(0, 1700)) +
        #annotate("text", x = 1, y = 1400, label = paste0("variance difference P: ", round(pval_df2$p.value, 4))) +
        theme_classic() +
        my_theme
    saveRDS(tissue_pcawg_vs_healthy_cnv_fig, file.path(plotdir, "pcawg", name, paste0("pcawg_", cancer_name, "_vs_healthy_cnv.rds")))
    ggsave(file.path(plotdir, "pcawg", name, paste0("pcawg_", cancer_name, "_vs_healthy_cnv.pdf")), tissue_pcawg_vs_healthy_cnv_fig)
    return(tissue_pcawg_vs_healthy_cnv)
}

create_pcawg_nuclear_sigs = function(pcawg_nuclear_sigs, high_mut_samples, name){
    
    
    # Create contri matrix for the correct samples.
    contri_m = pcawg_nuclear_sigs %>% dplyr::rename(sample_id = submitted_donor_id) %>% 
        dplyr::filter(sample_id %in% high_mut_samples) %>% 
        dplyr::select(-`Cancer Types`, -Accuracy) %>% 
        tibble::column_to_rownames("sample_id") %>% 
        t()
    
    # Remove absent signatures, so they don't use any colours
    contri_m = contri_m[rowSums(contri_m) != 0,]
    
    # Create plot
    pcawg_blood_nuclear_sig_fig = plot_contribution(contri_m)
    ggsave(file.path(plotdir, "pcawg", "high_mut_samples", paste0("pcawg_", name, "_nuclear_sig_fig.pdf")))
}


create_pcawg_drivers = function(pcawg_drivers, cancer_types, high_mut_samples, name){
    
    # Identify correct tumor types and which ones are high vs low mut.
    pcawg_drivers_mut = pcawg_drivers %>% 
        dplyr::filter(ttype %in% cancer_types) %>% 
        dplyr::mutate(mutation_count = ifelse(submitted_donor_id %in% high_mut_samples, "High", "Low"))
    
    # Look at number of drivers
    pcawg_drivers_ndrivers = pcawg_drivers_mut %>% 
        dplyr::group_by(submitted_donor_id) %>% 
        dplyr::summarise(n_drivers = dplyr::n(), ttype = ttype[[1]], mutation_count = mutation_count[[1]])
    high_mut_pcawg_driver_fig = ggplot(pcawg_drivers_ndrivers, aes(x = mutation_count, y = n_drivers, colour = ttype)) +
        geom_quasirandom() +
        labs(x = "Mitochondrial mutation count", y = "Nuclear drivers", colour = "Tumour type") +
        theme_BM() +
        my_theme
    ggsave(file.path(plotdir, "pcawg", "high_mut_samples", paste0("pcawg_", name, "_pcawg_drivers_amount.pdf")))
    
    
    # Create overview of which drivers are in which sample
    abundant_genes = pcawg_drivers_mut %>%
        dplyr::count(gene) %>% 
        dplyr::filter(n >= 10) %>% 
        pull(gene)
    pcawg_abundant_drivers_mut = pcawg_drivers_mut %>% 
        dplyr::filter(gene %in% abundant_genes)
    
    
    drivers_heat_fig = ggplot(pcawg_abundant_drivers_mut, aes(x = submitted_donor_id, y = gene)) +
        geom_raster() +
        facet_grid(. ~ mutation_count, scale = "free_x", space = "free_x") +
        theme_BM() +
        my_theme
    ggsave(file.path(plotdir, "pcawg", "high_mut_samples", paste0("pcawg_", name, "_pcawg_drivers_heat.pdf")))
    
    
    # Look at the mut ratio of driver genes.
    nr_samples = pcawg_drivers_mut %>%
        dplyr::filter(!duplicated(submitted_donor_id)) %>% 
        dplyr::count(mutation_count, name = "nr_samples")
    ratio_pcawg_drivers = pcawg_abundant_drivers_mut %>% 
        dplyr::group_by(gene, mutation_count, .drop = FALSE) %>% 
        dplyr::summarise(nr_mutated = dplyr::n(), .groups = "drop") %>%
        tidyr::complete(gene, mutation_count, fill = list(nr_mutated = 0)) %>%
        dplyr::left_join(nr_samples, by = "mutation_count") %>% 
        dplyr::mutate(ratio_mutated = nr_mutated / nr_samples)
    
    driver_ratio_fig = ggplot(ratio_pcawg_drivers, aes(x = gene, y = ratio_mutated, fill = mutation_count)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_BM() +
        my_theme
    ggsave(file.path(plotdir, "pcawg", "high_mut_samples", paste0("pcawg_", name, "_pcawg_drivers_ratio.pdf")))
}

