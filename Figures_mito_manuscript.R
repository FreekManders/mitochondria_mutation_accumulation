library(tidyverse)
library(gridExtra)
library(ggpubr)
 
# intestine_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "intestine_ageline.rds")) +
#     ggtitle("Intestine")
# colon_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "colon_ageline.rds")) +
#     ggtitle("Colon")
# blood_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "blood_ageline.rds")) + 
#     ggtitle("Blood")
# 
# coverage_fig = readRDS(file.path(plotdir, "coverage.rds"))

somatic_variant_df = read_tsv(file.path(r_wd, "somatic_variant.txt"))
annotated_variants = read_tsv(file.path(r_wd, "annotated_variants.txt")) # Filter out indels after this step.
total_sample_mutation_freq = read_tsv(file.path(r_wd, "total_sample_mutation_freq.txt"))

### Figure 1
cnv_df_all = readRDS(file.path(r_wd, "cnv.rds"))
blood_glm_m = readRDS(file.path(model_wd, "ageline.rds"))
colon_glm_m = readRDS(file.path(model_wd, "colon_ageline.rds"))
intestine_glm_m = readRDS(file.path(model_wd, "intestine_ageline.rds"))

coverage_fig = plot_coverage_per_patient(cnv_df_all)
blood_ageline_fig = plot_ageline(blood_glm_m) + 
    ggtitle("Blood") +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, 16.5), "points"))
colon_ageline_fig = plot_ageline(colon_glm_m) +
    ggtitle("Colon") +
    labs(y = "")
intestine_ageline_fig = plot_ageline(intestine_glm_m) +
    ggtitle("Intestine") +
    labs(y = "")



layout = matrix(c(1,1,1,2,3,4), nrow = 2, byrow = TRUE)
pdf(file.path("Manuscript_figures", "Figure1.pdf"), width = 6.8, height = 5, useDingbats = FALSE)
grid.arrange(coverage_fig, blood_ageline_fig, colon_ageline_fig, intestine_ageline_fig, layout_matrix = layout)
dev.off()


### Figure 2
ageline_cnv_model = readRDS(file.path(model_wd, "cnv_ageline.rds"))
colon_cnv_model = readRDS(file.path(model_wd, "colon_cnv_ageline.rds"))
intestine_cnv_model = readRDS(file.path(model_wd, "intestine_cnv_ageline.rds"))
blood_cnv_ageline_fig = plot_cnv_ageline(ageline_cnv_model) + 
    ggtitle("Blood")
colon_cnv_ageline_fig = plot_cnv_ageline(colon_cnv_model) + 
    ggtitle("Colon")
intestine_cnv_ageline_fig = plot_cnv_ageline(intestine_cnv_model) + 
    ggtitle("Intestine")

pdf(file.path("Manuscript_figures", "Figure2.pdf"), width = 6.8, height = 3, useDingbats = FALSE)
grid.arrange(blood_cnv_ageline_fig, colon_cnv_ageline_fig, intestine_cnv_ageline_fig, nrow = 1)
dev.off()

### Figure 3

blood_healthy_spectra_fig = readRDS("plot/pcawg/Healthy blood/tissue_pcawg_spectra/healthy_spectrum.rds") +
    theme(plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))
blood_cancer_spectra_fig = readRDS("plot/pcawg/Healthy blood/tissue_pcawg_spectra/Blood cancer_spectrum.rds") +
    theme(plot.margin = unit(c(0,5.5,0,5.5), "pt"))
colon_healthy_spectra_fig = readRDS("plot/pcawg/Normal colon/tissue_pcawg_spectra/healthy_spectrum.rds") +
    theme(plot.margin = unit(c(0,5.5,0,5.5), "pt"))
colon_cancer_spectra_fig = readRDS("plot/pcawg/Normal colon/tissue_pcawg_spectra/Colon cancer_spectrum.rds") +
    theme(plot.margin = unit(c(0,5.5,5.5,5.5), "pt"))

spectra_fig = ggarrange(blood_healthy_spectra_fig, blood_cancer_spectra_fig, 
          colon_healthy_spectra_fig, colon_cancer_spectra_fig, 
          common.legend = TRUE, nrow = 4, legend = "bottom", align = "v")

#blood_age_fig = readRDS("plot/pcawg/Normal blood/Blood cancerpcawg_age_model_fig.rds")
#colon_age_fig = readRDS("plot/pcawg/Normal colon/Colon cancerpcawg_age_model_fig.rds")

pcawg_blood_glm_m = readRDS("plot/pcawg/Healthy blood/Blood cancerpcawg_age_model.rds")
get_ci_params(pcawg_blood_glm_m)
pcawg_colon_glm_m = readRDS("plot/pcawg/Normal colon/Colon cancerpcawg_age_model.rds")
get_ci_params(pcawg_colon_glm_m)

blood_age_fig = plot_pcawg_model(pcawg_blood_glm_m)
colon_age_fig = plot_pcawg_model(pcawg_colon_glm_m)


layout = matrix(c(2,1,3,1), nrow = 2, byrow = TRUE)

pdf(file.path("Manuscript_figures", "Figure3.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
grid.arrange(spectra_fig, blood_age_fig, colon_age_fig, layout_matrix = layout)
dev.off()


### Figure 4
# dx1_cancer_m = readRDS(file.path(model_wd, "chemo_1st_cancer_ageline.rds"))
# get_ci_params(dx1_cancer_m)
# fu_m = readRDS(file.path(model_wd, "chemo_ageline_FU.rds"))
# get_ci_params(fu_m)
# ko_df = readRDS(file.path(r_wd,"ko.rds"))

cancer2_m = readRDS(file.path(model_wd, "chemo_2nd_cancer_ageline.rds"))
get_ci_params(cancer2_m)

dx1_fu_dx2_cancer_m = readRDS(file.path(model_wd, "chemo_dx1_fu_dx2_ageline.rds"))
get_ci_params(dx1_fu_dx2_cancer_m)

cb_chemo_freq = readRDS(file.path(r_wd,"CB_chemo.rds"))

cb_chemo_cnv = readRDS(file.path(r_wd, "CB_chemo_cnv.rds"))


cancer2_ageline_fig = plot_chemo_model(cancer2_m)
cb_chemo_fig = plot_cb_chemo(cb_chemo_freq) +
     theme(plot.margin = unit(c(5.5, 5.5, 5.5, 10), "points"))
dx1_fu_dx2_ageline_fig = plot_chemo_model(dx1_fu_dx2_cancer_m)
cnv_cb_chemo_fig = plot_cb_chemo_cnv(cb_chemo_cnv) + theme(plot.margin = unit(c(5.5,5.5,5.5,1), "pt"))

layout = matrix(c(1, 2, 3, 3, 4, 4), byrow = TRUE, nrow = 3)

pdf(file.path("Manuscript_figures", "Figure4.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
grid.arrange(dx1_fu_dx2_ageline_fig, cancer2_ageline_fig, 
             cb_chemo_fig, cnv_cb_chemo_fig, layout_matrix = layout)
dev.off()


### Figure 5

hsct_m = readRDS(file.path(model_wd, "hsct_ageline.rds"))
get_ci_params(hsct_m)
bloodtype_m = readRDS(file.path(model_wd, "cnv_bloodtype_transplant.rds"))

hsct_ageline_fig = plot_hsct_model(hsct_m) +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, 11.5), "points"))

bloodtype_hsct_cnv_fig = plot_cnv_model(bloodtype_m, var = "full_blood_type", col = "state", remove_guide = FALSE, plot_p = FALSE, size = 0.5) +
    labs(x = "HSPC location", color = "State") +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, 2.5), "points"))


layout = matrix(c(1,1,2,2), nrow = 2, byrow = TRUE)

pdf(file.path("Manuscript_figures", "Figure5.pdf"), width = 6.8, height = 4, useDingbats = FALSE)
grid.arrange(hsct_ageline_fig, bloodtype_hsct_cnv_fig, layout_matrix = layout)
dev.off()


### Figure S1
per_base_fig = readRDS(file.path(plotdir, "per_base_coverage", "mitochondrial_coverage_AC63.rds"))
#per_base_fig2 = readRDS(file.path(plotdir, "per_base_coverage", "mitochondrial_coverage_HSCT14.rds"))
somatic_variant_df = readRDS(file.path(r_wd, "somatic_variant_df.rds"))
indel_df = readRDS(file.path(r_wd, "indel_df.rds"))
cov_m = readRDS(file.path(model_wd, "sequencing_depth.rds"))
higher_lower_muts = readRDS(file.path(r_wd, "higher_lower_muts.rds"))
#liver_freq = readRDS(file.path(r_wd, "liver_freq.rds"))
freq_cnv_m = readRDS(file.path(model_wd, "freq_cnv.rds"))
somatic_variant_annotate_length = readRDS(file.path(r_wd, "somatic_variant_annotate_length.rds"))
somatic_variant_annotate = readRDS(file.path(r_wd, "somatic_variant_annotate.rds"))


vaf_fig = plot_vaf_position(somatic_variant_df)
indel_vaf_fig = plot_vaf_position(indel_df, muttype = "indel")
alt_fig = plot_alt_position(somatic_variant_df)
indel_alt_fig = plot_alt_position(indel_df, muttype = "indel")
snv_fig = ggarrange(vaf_fig, alt_fig, nrow = 2, common.legend = TRUE, legend = "bottom", align = "v")
indel_fig = ggarrange(indel_vaf_fig, indel_alt_fig, nrow = 2, common.legend = TRUE, legend = "bottom", align = "v")
gene_length_fig = plot_gene(somatic_variant_annotate_length, length_correction = TRUE) +
    theme(plot.margin = unit(c(5.5, 5.5, 0.0, 5.5), "points"))
effect_fig = plot_effect(somatic_variant_annotate) +
    theme(plot.margin = unit(c(0.0, 5.5, 8.5, 5.5), "points"))
seqdepth_fig = plot_cov_model2(cov_m) 
freq_cnv_fig = plot_freq_cnv_model(freq_cnv_m)
above_below_ageline_fig = plot_above_below_ageline(higher_lower_muts)
#liver_ageline_fig = plot_liver_ageline(liver_freq)

#agelines_sup_fig = ggarrange(seqdepth_fig, liver_ageline_fig, align = "v", common.legend = T, nrow = 2)
#agelines_sup_fig = plot_grid(seqdepth_fig, liver_ageline_fig, align = "v", axis = "rl", nrow = 2)

layout = matrix(c(1, 1, 2, 3, 2, 3, 4, 4, 5, 5, 6, 7, 8, 9), byrow = TRUE, nrow = 7)

pdf(file.path("Manuscript_figures", "Supplemental_figure1.pdf"), width = 6.8, height = 10, useDingbats = FALSE)
grid.arrange(per_base_fig, snv_fig, indel_fig, gene_length_fig, effect_fig, 
             seqdepth_fig, freq_cnv_fig, above_below_ageline_fig, 
             layout_matrix = layout, heights = c(0.3, 0.3, 0.3, 0.6, 0.4, 0.3, 0.3))
dev.off()

### Figure S2
cov_cnv_m = readRDS(file.path(model_wd, "cnv_cov.rds"))
#liver_cnv = readRDS(file.path(r_wd, "liver_cnv.rds"))


seqdepth_cnv_fig = plot_cnv_model(cov_cnv_m, var = "cov") + xlab("Sequencing depth")
#liver_cnv_ageline_fig = plot_cnv_liver_ageline(liver_cnv)

layout = matrix(c(1, 2), byrow = TRUE, nrow = 1)

pdf(file.path("Manuscript_figures", "Supplemental_figure2.pdf"), width = 6.8, height = 4, useDingbats = FALSE)
grid.arrange(seqdepth_cnv_fig, layout_matrix = layout)
dev.off()

### Figure S3
blood_cnv_fig = readRDS("plot/pcawg/Healthy blood/pcawg_Blood cancer_vs_healthy_cnv.rds")
colon_cnv_fig = readRDS("plot/pcawg/Normal colon/pcawg_Colon cancer_vs_healthy_cnv.rds")
mut_mat = readRDS(file.path(plotdir, "pcawg", "blood_colon_mut_mat.rds"))
blood_colon_profile_fig = plot_96_profile(mut_mat, condensed = TRUE) +
    spectrum_theme
aml_cnv = readRDS(file.path(r_wd, "aml_cnv.rds"))
aml_cnv_fig = plot_aml_cnv_model(aml_cnv)

layout = matrix(c(1, 1, 2, 3, 4, 5), byrow = TRUE, nrow = 3)
pdf(file.path("Manuscript_figures", "Supplemental_figure3.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
grid.arrange(blood_colon_profile_fig, blood_cnv_fig, colon_cnv_fig, aml_cnv_fig, layout_matrix = layout)
dev.off()

### Figure S4
chemo_cancer_cnv_m = readRDS(file.path(model_wd, "cnv_chemo_2ndcancer.rds"))
chemo_cnv_fig = plot_cnv_model(chemo_cancer_cnv_m, plot_p = FALSE)


#layout = matrix(c(1, 2, 3, 3, 3, 3), byrow = TRUE, nrow = 3)
pdf(file.path("Manuscript_figures", "Supplemental_figure4.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
chemo_cnv_fig
dev.off()


### Figure S5
hsct_cnv_m = readRDS(file.path(model_wd, "cnv_hsct.rds"))
recipient_m = readRDS(file.path(model_wd, "cnv_hsct_time.rds"))

mean_recipient_fig = plot_mean_hsct_cnv(hsct_cnv_m$data) + theme(plot.margin = unit(c(5.5,5.5,5.5,10), "pt"))
recipient_fig = plot_recipient_time_model(recipient_m)

pdf(file.path("Manuscript_figures", "Supplemental_figure5.pdf"), width = 6.8, height = 4, useDingbats = FALSE)
ggarrange(mean_recipient_fig, recipient_fig, nrow = 2)
dev.off()




# 
# library(tidyverse)
# library(gridExtra)
# library(ggpubr)
# 
# # intestine_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "intestine_ageline.rds")) +
# #     ggtitle("Intestine")
# # colon_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "colon_ageline.rds")) +
# #     ggtitle("Colon")
# # blood_ageline_fig = readRDS(file.path(plotdir, "mut_accumulation", "blood_ageline.rds")) + 
# #     ggtitle("Blood")
# # 
# # coverage_fig = readRDS(file.path(plotdir, "coverage.rds"))
# 
# somatic_variant_df = read_tsv(file.path(r_wd, "somatic_variant.txt"))
# annotated_variants = read_tsv(file.path(r_wd, "annotated_variants.txt")) # Filter out indels after this step.
# total_sample_mutation_freq = read_tsv(file.path(r_wd, "total_sample_mutation_freq.txt"))
# 
# ### Figure 1
# cnv_df = read_tsv(file.path(r_wd, "cnv.txt"))
# blood_glm_m = readRDS(file.path(model_wd, "ageline.rds"))
# colon_glm_m = readRDS(file.path(model_wd, "colon_ageline.rds"))
# intestine_glm_m = readRDS(file.path(model_wd, "intestine_ageline.rds"))
# coverage_fig = plot_coverage_per_patient(cnv_df)
# blood_ageline_fig = plot_ageline(blood_glm_m) + 
#     ggtitle("Blood") +
#     theme(plot.margin = unit(c(5.5, 5.5, 5.5, 16.5), "points"))
# colon_ageline_fig = plot_ageline(colon_glm_m) +
#     ggtitle("Colon") +
#     labs(y = "")
# intestine_ageline_fig = plot_ageline(intestine_glm_m) +
#     ggtitle("Intestine") +
#     labs(y = "")
# 
# 
# 
# layout = matrix(c(1,1,1,2,3,4), nrow = 2, byrow = TRUE)
# pdf(file.path("Manuscript_figures", "Figure1.pdf"), width = 6.8, height = 5, useDingbats = FALSE)
# grid.arrange(coverage_fig, blood_ageline_fig, colon_ageline_fig, intestine_ageline_fig, layout_matrix = layout)
# dev.off()
# 
# 
# ### Figure 2
# ageline_cnv_model = readRDS(file.path(model_wd, "cnv_ageline.rds"))
# colon_cnv_model = readRDS(file.path(model_wd, "colon_cnv_ageline.rds"))
# intestine_cnv_model = readRDS(file.path(model_wd, "intestine_cnv_ageline.rds"))
# blood_cnv_ageline_fig = plot_cnv_ageline(ageline_cnv_model) + 
#     ggtitle("Blood")
# colon_cnv_ageline_fig = plot_cnv_ageline(colon_cnv_model) + 
#     ggtitle("Colon")
# intestine_cnv_ageline_fig = plot_cnv_ageline(intestine_cnv_model) + 
#     ggtitle("Intestine")
# 
# pdf(file.path("Manuscript_figures", "Figure2.pdf"), width = 6.8, height = 3, useDingbats = FALSE)
# grid.arrange(blood_cnv_ageline_fig, colon_cnv_ageline_fig, intestine_cnv_ageline_fig, nrow = 1)
# dev.off()
# 
# ### Figure 3
# 
# blood_healthy_spectra_fig = readRDS("plot/pcawg/Healthy blood/tissue_pcawg_spectra/healthy_spectrum.rds") +
#     theme(plot.margin = unit(c(5.5,5.5,0,5.5), "pt"))
# blood_cancer_spectra_fig = readRDS("plot/pcawg/Healthy blood/tissue_pcawg_spectra/Blood cancer_spectrum.rds") +
#     theme(plot.margin = unit(c(0,5.5,0,5.5), "pt"))
# colon_healthy_spectra_fig = readRDS("plot/pcawg/Normal colon/tissue_pcawg_spectra/healthy_spectrum.rds") +
#     theme(plot.margin = unit(c(0,5.5,0,5.5), "pt"))
# colon_cancer_spectra_fig = readRDS("plot/pcawg/Normal colon/tissue_pcawg_spectra/Colon cancer_spectrum.rds") +
#     theme(plot.margin = unit(c(0,5.5,5.5,5.5), "pt"))
# 
# spectra_fig = ggarrange(blood_healthy_spectra_fig, blood_cancer_spectra_fig, 
#                         colon_healthy_spectra_fig, colon_cancer_spectra_fig, 
#                         common.legend = TRUE, nrow = 4, legend = "bottom", align = "v")
# 
# #blood_age_fig = readRDS("plot/pcawg/Normal blood/Blood cancerpcawg_age_model_fig.rds")
# #colon_age_fig = readRDS("plot/pcawg/Normal colon/Colon cancerpcawg_age_model_fig.rds")
# 
# pcawg_blood_glm_m = readRDS("plot/pcawg/Healthy blood/Blood cancerpcawg_age_model.rds")
# get_ci_params(pcawg_blood_glm_m)
# pcawg_colon_glm_m = readRDS("plot/pcawg/Normal colon/Colon cancerpcawg_age_model.rds")
# get_ci_params(pcawg_colon_glm_m)
# 
# blood_age_fig = plot_pcawg_model(pcawg_blood_glm_m)
# colon_age_fig = plot_pcawg_model(pcawg_colon_glm_m)
# 
# 
# layout = matrix(c(2,1,3,1), nrow = 2, byrow = TRUE)
# 
# pdf(file.path("Manuscript_figures", "Figure3.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
# grid.arrange(spectra_fig, blood_age_fig, colon_age_fig, layout_matrix = layout)
# dev.off()
# 
# 
# ### Figure 4
# dx1_cancer_m = readRDS(file.path(model_wd, "chemo_1st_cancer_ageline.rds"))
# get_ci_params(dx1_cancer_m)
# dx2_cancer_m = readRDS(file.path(model_wd, "chemo_2nd_cancer_ageline.rds"))
# get_ci_params(dx2_cancer_m)
# fu_m = readRDS(file.path(model_wd, "chemo_ageline_FU.rds"))
# get_ci_params(fu_m)
# cb_chemo_freq = readRDS(file.path(r_wd,"CB_chemo.rds"))
# hsct_m = readRDS(file.path(model_wd, "hsct_ageline.rds"))
# get_ci_params(hsct_m)
# ko_df = readRDS(file.path(r_wd,"ko.rds"))
# 
# 
# 
# dx1_cancer_ageline_fig = plot_chemo_model(dx1_cancer_m)
# dx2_cancer_ageline_fig = plot_chemo_model(dx2_cancer_m)
# fu_ageline_fig = plot_chemo_model(fu_m) +
#     theme(plot.margin = unit(c(5.5, 5.5, 5.5, 2.5), "points"))
# cb_chemo_fig = plot_cb_chemo(cb_chemo_freq)
# hsct_ageline_fig = plot_hsct_model(hsct_m)
# ko_fig = plot_ko(ko_df) +
#     theme(plot.margin = unit(c(5.5, 5.5, 5.5, 3.1), "points"))
# 
# pdf(file.path("Manuscript_figures", "Figure4.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
# grid.arrange(dx1_cancer_ageline_fig, fu_ageline_fig, 
#              dx2_cancer_ageline_fig, cb_chemo_fig, 
#              hsct_ageline_fig, ko_fig)
# dev.off()
# 
# 
# ### Figure 5
# chemo_cancer_cnv_m = readRDS(file.path(model_wd, "cnv_chemo_2ndcancer.rds"))
# hsct_cnv_m = readRDS(file.path(model_wd, "cnv_hsct.rds"))
# cb_chemo_cnv = readRDS(file.path(r_wd, "CB_chemo_cnv.rds"))
# #target_cnv = readRDS(file.path(r_wd, "target_cnv.rds"))
# cnv_df = readRDS(file.path(r_wd, "cnv_df_clean.rds"))
# 
# chemo_cnv_fig = plot_cnv_model(chemo_cancer_cnv_m, plot_p = FALSE)
# hsct_cnv_fig = plot_cnv_model(hsct_cnv_m)
# cb_chemo_fig = plot_cb_chemo_cnv(cb_chemo_cnv)
# #aml_bulk_mtreads_fig = compare_aml_bulk_cnv(target_cnv)
# ko_cnv_fig = plot_ko_cnv(cnv_df)
# 
# pdf(file.path("Manuscript_figures", "Figure5.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
# grid.arrange(chemo_cnv_fig, cb_chemo_fig, ko_cnv_fig, hsct_cnv_fig)
# dev.off()
# 
# 
# ### Figure S1
# per_base_fig = readRDS(file.path(plotdir, "per_base_coverage", "mitochondrial_coverage_AC63.rds"))
# #per_base_fig2 = readRDS(file.path(plotdir, "per_base_coverage", "mitochondrial_coverage_HSCT14.rds"))
# somatic_variant_df = readRDS(file.path(r_wd, "somatic_variant_df.rds"))
# indel_df = readRDS(file.path(r_wd, "indel_df.rds"))
# 
# vaf_fig = plot_vaf_position(somatic_variant_df)
# indel_vaf_fig = plot_vaf_position(indel_df, muttype = "indel")
# alt_fig = plot_alt_position(somatic_variant_df)
# indel_alt_fig = plot_alt_position(indel_df, muttype = "indel")
# snv_fig = ggarrange(vaf_fig, alt_fig, nrow = 2, common.legend = TRUE, legend = "bottom", align = "v")
# indel_fig = ggarrange(indel_vaf_fig, indel_alt_fig, nrow = 2, common.legend = TRUE, legend = "bottom", align = "v")
# 
# layout = matrix(c(1, 1, 2, 3, 2, 3), byrow = TRUE, nrow = 3)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure1.pdf"), width = 6.8, height = 6, useDingbats = FALSE)
# grid.arrange(per_base_fig, snv_fig, indel_fig, layout_matrix = layout)
# dev.off()
# 
# ### Figure S2
# cov_m = readRDS(file.path(model_wd, "sequencing_depth.rds"))
# cov_cnv_m = readRDS(file.path(model_wd, "cnv_cov.rds"))
# 
# seqdepth_fig = plot_cov_model2(cov_m)
# seqdepth_cnv_fig = plot_cnv_model(cov_cnv_m, var = "cov") + xlab("Sequencing depth")
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure2.pdf"), width = 6.8, height = 2, useDingbats = FALSE)
# grid.arrange(seqdepth_fig, seqdepth_cnv_fig, nrow = 1)
# dev.off()
# 
# ### Figure S3
# higher_lower_muts = readRDS(file.path(r_wd, "higher_lower_muts.rds"))
# above_below_ageline_fig = plot_above_below_ageline(higher_lower_muts)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure3.pdf"), width = 3.4, height = 2, useDingbats = FALSE)
# above_below_ageline_fig
# dev.off()
# 
# ### Figure S4
# liver_freq = readRDS(file.path(r_wd, "liver_freq.rds"))
# liver_cnv = readRDS(file.path(r_wd, "liver_cnv.rds"))
# 
# liver_ageline_fig = plot_liver_ageline(liver_freq)
# liver_cnv_ageline_fig = plot_cnv_liver_ageline(liver_cnv)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure4.pdf"), width = 6.8, height = 4, useDingbats = FALSE)
# ggarrange(liver_ageline_fig, liver_cnv_ageline_fig, nrow = 2, align = "v")
# dev.off()
# 
# ### Figure S5
# freq_cnv_m = readRDS(file.path(model_wd, "freq_cnv.rds"))
# freq_cnv_fig = plot_freq_cnv_model(freq_cnv_m)
# pdf(file.path("Manuscript_figures", "Supplemental_figure5.pdf"), width = 3.4, height = 2, useDingbats = FALSE)
# freq_cnv_fig
# dev.off()
# 
# 
# 
# ### Figure S6
# blood_cnv_fig = readRDS("plot/pcawg/Normal blood/pcawg_Blood cancer_vs_healthy_cnv.rds")
# colon_cnv_fig = readRDS("plot/pcawg/Normal colon/pcawg_Colon cancer_vs_healthy_cnv.rds")
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure6.pdf"), width = 6.8, height = 2, useDingbats = FALSE)
# ggarrange(blood_cnv_fig, colon_cnv_fig, nrow = 1)
# dev.off()
# 
# 
# ### Figure S7
# foetus_freq = readRDS(file.path(r_wd, "foetus_freq.rds"))
# fetus_cnv_m = readRDS(file.path(model_wd, "cnv_fetus_divstri.rds"))
# 
# fetus_ageline_fig2 = plot_fetus_model2(foetus_freq)
# fetus_cnv_fig = plot_cnv_model(fetus_cnv_m, plot_p = FALSE)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure7.pdf"), width = 6.8, height = 2, useDingbats = FALSE)
# ggarrange(fetus_ageline_fig2, fetus_cnv_fig, nrow = 1)
# dev.off()
# 
# ### Figure S8
# somatic_variant_df = readRDS(file.path(r_wd, "somatic_variant_df.rds"))
# vaf_age_fig = plot_vaf_age_cor(somatic_variant_df)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure8.pdf"), width = 3.4, height = 2, useDingbats = FALSE)
# vaf_age_fig
# dev.off()
# 
# ### Figure S9
# somatic_variant_annotate_length = readRDS(file.path(r_wd, "somatic_variant_annotate_length.rds"))
# somatic_variant_annotate = readRDS(file.path(r_wd, "somatic_variant_annotate.rds"))
# 
# gene_length_fig = plot_gene(somatic_variant_annotate_length, length_correction = TRUE)
# effect_fig = plot_effect(somatic_variant_annotate)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure9.pdf"), width = 6.8, height = 5, useDingbats = FALSE)
# ggarrange(gene_length_fig, effect_fig, nrow = 2, heights = c(0.6, 0.4))
# dev.off()
# 
# 
# ### Figure S10
# aml_cnv = readRDS(file.path(r_wd, "aml_cnv.rds"))
# aml_cnv_fig = plot_aml_cnv_model(aml_cnv)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure10.pdf"), width = 3.4, height = 2, useDingbats = FALSE)
# aml_cnv_fig
# dev.off()
# 
# 
# ### Figure S11
# hsct_cnv_m = readRDS(file.path(model_wd, "cnv_hsct.rds"))
# recipient_m = readRDS(file.path(model_wd, "cnv_hsct_time.rds"))
# 
# mean_recipient_fig = plot_mean_hsct_cnv(hsct_cnv_m$data)
# recipient_fig = plot_recipient_time_model(recipient_m)
# 
# pdf(file.path("Manuscript_figures", "Supplemental_figure11.pdf"), width = 6.8, height = 4, useDingbats = FALSE)
# ggarrange(mean_recipient_fig, recipient_fig, nrow = 2)
# dev.off()







