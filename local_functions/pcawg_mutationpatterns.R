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
    pcawg_snv$effect_type = pcawg_snv$var_type
    
    plot_save_pcawg_spectra(pcawg_snv, variants, "pcawg", name, type = "all_pcawg", name)
    
    blood_pcawg_snv = pcawg_snv %>% 
        dplyr::filter(cancer_type %in% cancer_types)
    plot_save_pcawg_spectra(blood_pcawg_snv, variants, cancer_name, name, type = "tissue_pcawg", name)
    
    high_mut_snv = blood_pcawg_snv %>% 
        dplyr::filter(sample_id %in% high_mut_samples)
    plot_save_pcawg_spectra(high_mut_snv, variants, paste0(cancer_name, "_high_mut"), name, type = "high_mut_pcwag", name)
    
    low_mut_snv = blood_pcawg_snv %>% 
        dplyr::filter(!sample_id %in% high_mut_samples)
    plot_save_pcawg_spectra(high_mut_snv, low_mut_snv, paste0(cancer_name, ": High"), paste0(cancer_name, ": Low"), type = "high_vs_low_mut_pcawg", name)
}

plot_save_pcawg_spectra = function(pcawg_variants, healthy_variants, name1, name2, type, dir_name){
    
    if (type == "high_vs_low_mut_pcawg"){
        remove_duplis1 = FALSE
    } else{
        remove_duplis1 = TRUE
        healthy_variants = dplyr::rename(healthy_variants, sample_id = sample)
    }
    
    out_dir = file.path(plotdir, "pcawg", dir_name, paste0(type, "_spectra"))
    if (!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    old_dir = setwd(out_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Number of overlapping mutations between the two groups.
    unique_pcawg = unique(pcawg_variants[,c("mut_pos", "ref", "alt")])
    unique_healhty = unique(healthy_variants[,c("mut_pos", "ref", "alt")])
    overlap_tbl = inner_join(unique_pcawg, unique_healhty, by = c("mut_pos", "ref", "alt"))
    write_tsv(overlap_tbl, file.path(r_wd, paste0(name1, type, "_overlapping_mutations.txt")))
    
    # pcawg spectrum
    pcawg_spectrum_fig = plot_mt_spectrum(pcawg_variants, ref_genome, name1, remove_duplis = FALSE) +
        spectrum_theme
    ggsave(paste0(name1, "_spectrum.pdf"), pcawg_spectrum_fig)
    saveRDS(pcawg_spectrum_fig, paste0(name1, "_spectrum.rds"))
    
    # pcawg compare spectra
    mito_cossims = mito_cossim_spectra(healthy_variants, pcawg_variants, 
                                       remove_duplis1 = remove_duplis1, remove_duplis2 = FALSE)
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
    healthy_gr = to_granges_mt(healthy_variants, remove_duplis = remove_duplis1)
    pcawg_gr = to_granges_mt(pcawg_variants, remove_duplis = FALSE)
    
    grl = GRangesList(healthy_gr, pcawg_gr)
    names(grl) = c(name2, name1)
    mut_mat = mut_matrix(grl, ref_genome)
    pcawg_profile_comparison_fig = plot_compare_profiles(mut_mat[,1], mut_mat[,2], profile_names = names(grl))
    ggsave(paste0(name1, "_profile_comparison.pdf"), pcawg_profile_comparison_fig)
    
    
    signatures = get_known_signatures()
    #context_i = str_detect(rownames(mut_mat), "C>T|T>C")
    #mut_mat = mut_mat[context_i,]
    #signatures = signatures[context_i,]
    
    fit_res = fit_to_signatures(mut_mat, signatures)
    plot_contribution(fit_res$contribution)
    strict_refit = fit_to_signatures_strict(mut_mat, signatures)
    plot_contribution(strict_refit$fit_res$contribution)
    boot_refit = fit_to_signatures_bootstrapped(mut_mat, signatures)
    plot_bootstrapped_contribution(boot_refit$contribution, plot_type = "dotplot")
    
    cos_sim_m = cos_sim_matrix(mut_mat, signatures)
    plot_cosine_heatmap(cos_sim_m)
    
    # Determine distribution across mitochondrial genome
    pcawg_pos = pcawg_variants %>% 
        dplyr::mutate(group = name1) %>% 
        dplyr::select(group, sample_id, mut_pos)
    healthy_pos = healthy_variants %>% 
        dplyr::mutate(group = name2) %>% 
        dplyr::select(group, sample_id, mut_pos)
    if (remove_duplis1){
        healthy_pos = unique(healthy_pos)
    }
    pos_tbl = rbind(healthy_pos, pcawg_pos) %>% 
        dplyr::mutate(group = factor(group, levels = unique(group)))
    saveRDS(pos_tbl, paste0(name1, "_position_tbl.rds"))
    
    distribution_fig = plot_bin_distribution(pos_tbl)
    ggsave(paste0(name1, "_position.pdf"), distribution_fig)
    
    # Compare mutation effects
    pcawg_effects = pcawg_variants %>% 
        dplyr::mutate(
            effect_type = case_when(
                effect_type %in% c("nsSNP", "missense_variant", "start_lost", "stop_gained", "stop_lost") ~ "Missense",
                effect_type %in% c("synSNP", "synonymous_variant", "stop_retained_variant") ~ "Synonymous",
                effect_type %in% c("intergenic", "LocInfo", "ncUTR") ~ "Non-coding",
                is.na(var_type) ~ "Non-coding"
            ), 
            effect_type = factor(effect_type, levels = c("Missense", "Synonymous", "Non-coding")),
            group = name1) %>% 
        dplyr::select(group, effect_type)
    
    healthy_effects = healthy_variants %>% 
        dplyr::mutate(
            effect_type = case_when(
                effect_type %in% c("nsSNP", "missense_variant", "start_lost", "stop_gained", "stop_lost") ~ "Missense",
                effect_type %in% c("synSNP", "synonymous_variant", "stop_retained_variant") ~ "Synonymous",
                effect_type %in% c("intergenic", "LocInfo", "ncUTR") ~ "Non-coding",
                is.na(effect_type) ~ "Non-coding"
            ), 
            effect_type = factor(effect_type, levels = c("Missense", "Synonymous", "Non-coding")),
            group = name2) %>% 
        dplyr::select(group, effect_type)
    
    effects_tbl = rbind(healthy_effects, pcawg_effects) %>% 
        dplyr::mutate(effect_type = ifelse(effect_type == "Missense", "Missense", "Other"),
                      group = factor(group, levels = unique(group))) %>% 
        dplyr::group_by(group, effect_type, .drop = FALSE) %>%
        dplyr::count() %>% 
        dplyr::ungroup()
    saveRDS(effects_tbl, paste0(name1, "_effect_type_tbl.rds"))
    
    effect_m = effects_tbl %>% 
        pivot_wider(names_from = effect_type, values_from = n) %>%
        as.data.frame() %>% 
        column_to_rownames("group") %>% 
        as.matrix()

    effect_type_fig = plot_effect_type(effects_tbl)
    ggsave(paste0(name1, "_effect_type_comparison.pdf"), effect_type_fig)
    
    set.seed(1)
    chisq_res = effect_m %>% 
        chisq.test(simulate.p.value = TRUE) %>% 
        broom::tidy()
    write_tsv(chisq_res, file.path(r_wd, paste0(name1, type, "_effect_type_test.txt")))
    
    invisible(0)
}

plot_bin_distribution = function(pos_tbl){
    fig = ggplot(pos_tbl, aes(x = mut_pos, fill = group)) +
        geom_histogram(binwidth = 1000) +
        coord_cartesian(xlim = c(1, 16569)) +
        facet_grid(group ~ .) +
        scale_fill_manual(values = c("Healthy blood" = "#E64B35FF",
                                     "Normal colon" = "#E64B35FF",
                                     "Blood cancer" = "darkred",
                                     "Colon cancer" = "darkred",
                                     "Blood cancer: High" = "darkred",
                                     "Blood cancer: Low" = "darkred",
                                     "Colon cancer: High" = "darkred",
                                     "Colon cancer: Low" = "darkred"),
                          limits = force) +
        guides(fill = "none") +
        labs(x = "mtDNA position (bp)", y = "Base substitutions", fill = "Type") +
        theme_BM() +
        my_theme
    return(fig)
    # ggplot(pos_tbl, aes(x = mut_pos, color = group)) +
    #     geom_density() +
    #     coord_cartesian(xlim = c(1, 16569)) +
    #     theme_BM() +
    #     my_theme
}

plot_effect_type = function(effects_tbl){
    effect_m = effects_tbl %>% 
        pivot_wider(names_from = effect_type, values_from = n) %>%
        as.data.frame() %>% 
        column_to_rownames("group") %>% 
        as.matrix()
    
    set.seed(1)
    chisq_res = effect_m %>% 
        chisq.test(simulate.p.value = TRUE) %>% 
        broom::tidy()
    
    fig = ggplot(effects_tbl, aes(x = group, y = n, fill = effect_type)) +
        geom_bar(stat = "identity") +
        labs(x = "", y = "Base substitutions", fill = "Effect") +
        annotate("text", x = 1.5, y = 1.05*max(rowSums(effect_m)), label = paste0("P: ", round(chisq_res$p.value, 3)), size = 2) +
        scale_fill_manual(values = c("#fed9b7", "#00afb9")) +
        theme_BM() +
        my_theme
    return(fig)
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
    high_mut_samples = pcawg_age_freq$sample[pcawg_age_freq$freq >= 4]
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
    #histo_high_low_m = histo_high_low_m[,colSums(histo_high_low_m) >= 10]
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
    
    
    #pval_df = broom::tidy(wilcox.test(copy_number ~ type, data = blood_pcawg_vs_healthy_cnv))
    pval_df2 = broom::tidy(fligner.test(copy_number ~ type, data = tissue_pcawg_vs_healthy_cnv))
    
    tissue_pcawg_vs_healthy_cnv = tissue_pcawg_vs_healthy_cnv %>% 
        dplyr::mutate(type = factor(type, levels = c(name, cancer_name)))
    saveRDS(tissue_pcawg_vs_healthy_cnv, file.path(plotdir, "pcawg", name, paste0("pcawg_", cancer_name, "_vs_healthy_cnv.rds")))
    
    tissue_pcawg_vs_healthy_cnv_fig = plot_pcawg_healthy_vs_cnv(tissue_pcawg_vs_healthy_cnv)
    ggsave(file.path(plotdir, "pcawg", name, paste0("pcawg_", cancer_name, "_vs_healthy_cnv.pdf")), tissue_pcawg_vs_healthy_cnv_fig)
    return(tissue_pcawg_vs_healthy_cnv)
}

plot_pcawg_healthy_vs_cnv = function(tissue_pcawg_vs_healthy_cnv){
    m = lme(copy_number ~ type, random =  ~1 | patient, data = tissue_pcawg_vs_healthy_cnv)
    pval = summary(m)$tTable[2, "p-value"]
    pval <- formatC(pval,
                    digits = 4,
                    format = 'f')
    
    fig = ggplot(tissue_pcawg_vs_healthy_cnv, aes(x = type, y = copy_number, fill = type)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.25) +
        geom_quasirandom(groupOnX = TRUE, size = 0.5) +
        annotate("text", x = 1.5, y = 1700, label = paste0("P: ", pval), size = 2) +
        labs(fill = "Type", x = "", y = "mtDNA Copy Number") +
        coord_cartesian(ylim = c(0, 1700)) +
        scale_fill_manual(values = c("Healthy blood" = "#E64B35FF",
                                     "Normal colon" = "#E64B35FF",
                                     "Normal liver" = "#E64B35FF",
                                     "Blood cancer" = "darkred",
                                     "Colon cancer" = "darkred",
                                     "Liver cancer" = "darkred"), limits = force) +
        guides(fill = "none") +
        #annotate("text", x = 1, y = 1400, label = paste0("variance difference P: ", round(pval_df2$p.value, 4))) +
        theme_classic() +
        my_theme
    return(fig)
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

