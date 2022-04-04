library(tidyverse)
library(ggbeeswarm)
library(MutationalPatterns)
library(VariantAnnotation)
library(ggpubr)
library(ggrepel)
library(lme4)
library(nlme)
library(regioneR)
library(karyoploteR)
library(magrittr)
library(car)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gridExtra)
library(RColorBrewer)
library(devtools)
library(biomaRt)
library(extrafont)
library(ggeffects)
library(DESeq2)
library(biomaRt)
load_all("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/MutationalPatterns/")


# copy_data_from_hpc(source_dir = "~/hpc/pmc_vanboxtel/_old_structure/projects/Freek_mito/",
#                    dest_dir = "/Users/freekmanders/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/mito",
#                    sample_folder = "AML/",
#                    sample_reference_df = sample_reference_df)

#Questions: glm with identity link?
# Why does glm with mixed not work?
# Model selection and fitness tests of mito ageline.
# Why do BIC and AIC sometimes say the models are trained on different data.
# model selection fetus cnv.

#Source functions
function_wd <- "~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/mito/scripts/local_functions/"
functions_list <- list.files(function_wd, pattern = '.R$', ignore.case = T, full.names = T)
purrr::walk(functions_list, source)



# Ensure plot aesthetics look good.
# remotes::install_version("Rttf2pt1", version = "1.3.8") # Run only first time
# font_import() # Run only first time
loadfonts(device = "pdf")

my_theme = theme(text = element_text(size = 6, family = "Arial"),
                 legend.background = element_rect(fill="transparent", colour=NA),
                 legend.key = element_rect(fill="transparent", colour=NA),
                 axis.text = element_text(size = rel(1), 
                                          colour = "black"),
                 axis.ticks.x = element_line(colour = "black", size = 0.25),
                 axis.ticks.y = element_line(colour = "black", size = 0.25),
                 axis.line = element_line(size = 0.25),
                 strip.background = element_blank(),
                 plot.title = element_text(hjust = 0.5),
                 axis.ticks.length=unit(0.05, "cm"),
                 legend.key.size = unit(0.3, 'cm'),
                 plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
                 legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)

spectrum_theme = theme(text = element_text(size = 6, family = "Arial"),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      axis.text = element_text(size = rel(1), 
                               colour = "black"),
      axis.line = element_line(size = 0.25),
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(size = 0.5),
      panel.border = element_rect(size = 0.25),
      legend.key.size = unit(0.3, 'cm'),
      plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)
#Basic settings
main_dir = "~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Freek/mito"
plotdir <- file.path(main_dir, "plot")
vcf_wd <- file.path(main_dir, "annotated_vcfs")
wgs_wd <- file.path(main_dir, "output")
r_wd <- file.path(main_dir, "R_output")
model_wd <- file.path(main_dir, "models")
vcf_filter <- "*\\_snpSift.vcf$|*hg38.sorted.vcf$" # To select annotated VCFs
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
low_vaf_filter <- 0.01 # The lower and upper bound of VAF filtering, as by GATK recommendation
high_vaf_filter <- 1
setwd(main_dir)
sample_reference_df = read_tsv("reference_mito_table.txt", comment = "#") %>% 
  dplyr::mutate(state_name = case_when(state == "healthy" ~ "Healthy blood",
                                       state == "healthy_colon" ~ "Normal colon",
                                       state == "healthy_intestine" ~ "Healthy intestine",
                                       state == "healthy_liver" ~ "Normal liver",
                                       state == "KO" ~ "DNA repair knockout",
                                       state == "recipient" ~ "HSCT recipient",
                                       state == "trisomy" ~ "Trisomy 21",
                                       state == "FU" ~ "Follow-up",
                                       state == "DX1" ~ "Diagnosis",
                                       state == "DX2" ~ "Diagnosis 2",
                                       state == "CB_Chemo" ~ "In-vitro Chemotherapy",
                                       (state == "AML" | state == "ALL") & chemo == "DX1_cancer" ~ "Leukemia",
                                       (state == "AML" | state == "ALL") & chemo == "DX2_cancer" ~ "2nd Leukemia",
                                       str_detect(state, "bulk") ~ "Bulk control",
                                       TRUE ~ state),
                state_name = factor(state_name, levels = c("Healthy blood", "Normal colon", "Healthy intestine", "Normal liver", "Diagnosis",
                                    "Leukemia", "Follow-up", "Diagnosis 2", "2nd Leukemia", "In-vitro Chemotherapy", "HSCT recipient", "Bulk control")))

# Remove TARGET data that had bad nuclear results
bad_target_donors = c("PANVPB", "PARHRS", "PARLSL", "PARNAW", "PARZIA", "PASFLG", "PASIGA",
                      "PASNKZ", "PASTZK", "PASYEJ", "PATAIJ", "PATHIU", "PATJMY", "PATKBK",
                      "PATKWH")
sample_reference_df <- dplyr::filter(sample_reference_df, !patient %in% bad_target_donors)

# Calculate and plot cnv
cnv_df <- calculate_mito_cnv(wd = wgs_wd, 
                             sample_reference_df = sample_reference_df)
write_tsv(cnv_df, file.path(r_wd, "cnv.txt"))
saveRDS(cnv_df, file.path(r_wd, "cnv.rds"))
coverage_fig = plot_coverage_per_patient(cnv_df)
ggsave(file.path(plotdir, "coverage.pdf"), coverage_fig, width = 20)
saveRDS(coverage_fig, file.path(plotdir, "coverage.rds"))

#Per base coverage figures per patient
per_base_coverage <- calculate_mt_coverage(wd = wgs_wd, 
                                           mito_reference_df = sample_reference_df)
# Coverage figures
if(!dir.exists("plot/per_base_coverage")){
  dir.create("plot/per_base_coverage")
}
fig_l = lapply(unique(sample_reference_df$patient), function(x) { 
  plot_cov(per_base_coverage, x, "plot/per_base_coverage/")
})

mean(cnv_df$mt_mean)


# Remove low cov samples
removed_samples <- cnv_df[cnv_df$mt_mean < 1000,]$sample
write_tsv(as.data.frame(removed_samples), file.path(r_wd, "removed_samples.txt"))
sample_reference_df <- dplyr::filter(sample_reference_df, !(sample %in% removed_samples))
cnv_df <- dplyr::filter(cnv_df, !(sample %in% removed_samples))
per_base_coverage <- dplyr::filter(per_base_coverage, !(sample %in% removed_samples))


# Remove patients that have only one sample
nr_samples_per_patient = table(sample_reference_df$patient)
single_sample_patients = names(nr_samples_per_patient[nr_samples_per_patient == 1])
removed_samples2 = dplyr::filter(sample_reference_df, patient %in% single_sample_patients) %>% 
  pull(sample)
write_tsv(as.data.frame(removed_samples2), file.path(r_wd, "removed_samples_single_patient.txt"))
sample_reference_df <- dplyr::filter(sample_reference_df, !(patient %in% single_sample_patients))
cnv_df <- dplyr::filter(cnv_df, !(patient %in% single_sample_patients))
saveRDS(cnv_df, file.path(r_wd, "cnv_df_clean.rds"))
per_base_coverage <- dplyr::filter(per_base_coverage, !(patient %in% single_sample_patients))

#Correlation between read coverage of samples
per_base_cov_m = per_base_coverage %>% 
  dplyr::select(sample, coverage, position) %>% 
  pivot_wider(names_from = sample, values_from = coverage) %>% 
  dplyr::select(-position) %>% 
  as.matrix()
per_base_cossim = cos_sim_matrix(per_base_cov_m, per_base_cov_m)
median(per_base_cossim)
range(per_base_cossim)
per_base_cor = cor(per_base_cov_m)

# Calculate average read depth of mt.
mean_depth = mean(cnv_df$mt_mean)
n = nrow(cnv_df)
error = qt(0.975,df=n-1)*sd(cnv_df$mt_mean)/sqrt(n)
quantiles = quantile(cnv_df$mt_mean)
tibble("mean_depth" = mean_depth,
       "lower" = mean_depth - error,
       "upper" = mean_depth + error,
       "lower_25" = quantiles[2],
       "upper_25" = quantiles[4])

# Get somatic variants
somatic_variant_df = get_somatic_variant_df(ref_genome, vcf_wd, vcf_filter, sample_reference_df, low_vaf_filter, high_vaf_filter) %>% 
  dplyr::filter(!variant %in% c("MT:1700_T/C", "MT:567_A/AC", "MT:3397_A/G", "MT:7961_T/C", "COSM7419795")) # Remove several variants that are FP


# Filter out AHH1 variants that did not occur in the subclones.
# not_ahh1 = dplyr::filter(somatic_variant_df, patient != "AHH1")
# subclones = c("AHH1WTG2SCG6", "AHH1HPRT1CG6SCB8", "AHH1XPC3G2SC1G4", "AHH1MSH23H5SCE3", "AHH1UNG4F2SCG9")
# ahh1_subclones = dplyr::filter(somatic_variant_df, patient == "AHH1" & sample %in% subclones & overall_freq == 1)
# somatic_variant_df = rbind(not_ahh1, ahh1_subclones)
write_tsv(somatic_variant_df, file.path(r_wd, "somatic_variant.txt"))

# Get nr of missense, synonymous, ect.
annotated_variants <- grab_variant_annotations(vcf_filter = vcf_filter)
write_tsv(annotated_variants, file.path(r_wd, "annotated_variants.txt"))

# Get indel variants
indel_df <- filter(somatic_variant_df, 
                   nchar(somatic_variant_df$ref) > 1 | nchar(somatic_variant_df$alt) > 1)
indel_df$indel <- ifelse(nchar(indel_df$ref) > nchar(indel_df$alt), 
                         'deletion', 
                         'insertion')
indel_annotate <- left_join(indel_df, annotated_variants,
                            by = c('variant', 'sample'))


# Remove INDELs from somatic_variant_df
somatic_variant_df <- anti_join(somatic_variant_df, indel_df, 
                                by = 'variant')
somatic_variant_annotate <- left_join(somatic_variant_df, annotated_variants, 
                                      by = c('variant', 'sample'))

saveRDS(somatic_variant_df, file.path(r_wd, "somatic_variant_df.rds"))
saveRDS(indel_df, file.path(r_wd, "indel_df.rds"))
saveRDS(somatic_variant_annotate, file.path(r_wd, "somatic_variant_annotate.rds"))

#somatic_variant_df=readRDS(file.path(r_wd, "somatic_variant_df.rds"))

# Calculate number of mutations per sample
total_sample_mutation_freq <- count_variant_occurences(combined_mito_df = somatic_variant_df,
                                                       mito_reference_df = sample_reference_df) %>% 
  dplyr::mutate(state = factor(state, levels = unique(state)),
                bulk = factor(bulk, levels = unique(bulk)),
                cov = factor(cov, levels = unique(cov)),
                KO = factor(KO, levels = unique(KO)),
                source = factor(source, levels = unique(source)))

write_tsv(total_sample_mutation_freq, file.path(r_wd, "total_sample_mutation_freq.txt"))

# Create overview of mutation numbers
nr_muts_overview_fig = plot_nr_muts_overview(total_sample_mutation_freq)
ggsave(file.path(plotdir, "nr_muts_overview.pdf"), nr_muts_overview_fig, width = 12, height = 12)

nr_indels_overview_fig = plot_nr_indels_overview(indel_df)
ggsave(file.path(plotdir, "nr_indels_overview.pdf"), nr_indels_overview_fig)

# Subset DFs to only foetal variants

#Create spectra and profiles
healthy_variants <- subset(somatic_variant_annotate, 
                           state == 'healthy' &
                             !(patient == 'AHH1') &
                             bulk == 'clone')
plot_save_spectra(healthy_variants, "healthy")

colon_variants <- subset(somatic_variant_annotate, 
                           state == 'healthy_colon' &
                             !(patient == 'AHH1') &
                             bulk == 'clone')
plot_save_spectra(healthy_variants, "Colon")

intestine_variants <- subset(somatic_variant_annotate, 
                           state == 'healthy_intestine' &
                             !(patient == 'AHH1') &
                             bulk == 'clone')
plot_save_spectra(healthy_variants, "Small intestine")

liver_variants <- subset(somatic_variant_annotate, 
                           state == 'healthy_liver' &
                             !(patient == 'AHH1') &
                             bulk == 'clone')
plot_save_spectra(healthy_variants, "Liver")

aml_variants <- subset(somatic_variant_annotate, 
                       state == 'AML' &
                         bulk == 'clone' & 
                         (chemo == "DX1_cancer" | chemo == "No"))
plot_save_spectra(aml_variants, "aml")

hsct_recipient_variants <- subset(somatic_variant_annotate, 
                                  hsct == "recipient")
plot_save_spectra(hsct_recipient_variants, "hsct_recipient")

dx_variants = dplyr::filter(somatic_variant_annotate,
                                  chemo == "DX1")
plot_save_spectra(dx_variants, "DX1")

dx2_variants = dplyr::filter(somatic_variant_annotate,
                            chemo == "DX2")
plot_save_spectra(dx2_variants, "DX2")

cancer_variants = dplyr::filter(somatic_variant_annotate,
                                 chemo == "DX1_cancer")
plot_save_spectra(cancer_variants, "1st_cancer")

cancer2_variants = dplyr::filter(somatic_variant_annotate,
                             chemo == "DX2_cancer")
plot_save_spectra(cancer2_variants, "2nd_cancer")

chemo_fu_variants = dplyr::filter(somatic_variant_annotate,
                                  chemo == "FU")
plot_save_spectra(chemo_fu_variants, "FU")

chemo_cb_variants = dplyr::filter(somatic_variant_annotate,
                                  state == "CB_Chemo")
plot_save_spectra(chemo_cb_variants, "chemo_CB")


# Plot indel spectrum
indel_fig = plot_mt_indel_spectrum(indel_df, ref_genome)
ggsave(file.path(plotdir, "spectrum_indels.pdf"), indel_fig)

# Plot correlation vaf and age
vaf_age_cor_fig = plot_vaf_age_cor(healthy_variants)
ggsave(file.path(plotdir, "vaf_age_cor.pdf"), vaf_age_cor_fig)


#Rainfall plot
plot_save_rainfall(healthy_variants, "healthy")
plot_save_rainfall(colon_variants, "Colon")
plot_save_rainfall(intestine_variants, "Small intestine")
plot_save_rainfall(liver_variants, "Liver")
plot_save_rainfall(aml_variants, "aml")
plot_save_rainfall(hsct_recipient_variants, "hsct")
plot_save_rainfall(trisomy_variants, "trisomy")
plot_save_rainfall(dx_variants, "dx1")
plot_save_rainfall(dx2_variants, "dx2")
plot_save_rainfall(cancer_variants, "1st_cancer")
plot_save_rainfall(cancer2_variants, "2nd_cancer")
plot_save_rainfall(chemo_fu_variants, "chemo_fu")
plot_save_rainfall(chemo_cb_variants, "chemo_CB")


# Plot vaf and position
quantile(somatic_variant_df$VAF)
sd(somatic_variant_df$VAF)
quantile(somatic_variant_df$alt_ad)
sd(somatic_variant_df$alt_ad)
quantile(indel_df$VAF)
sd(indel_df$VAF)
quantile(indel_df$alt_ad)
sd(indel_df$alt_ad)


vaf_fig = plot_vaf_position(somatic_variant_df)
ggsave(file.path(plotdir, "vaf_position.pdf"), vaf_fig)

# Plot alt and position
alt_fig = plot_alt_position(somatic_variant_df)
ggsave(file.path(plotdir, "alt_position.pdf"), alt_fig)

# Plot snpeff annotation effect (HIGH/MODERATE/LOW/NONE)
effect_fig = plot_effect(somatic_variant_annotate)
ggsave(file.path(plotdir, "effect.pdf"), effect_fig, width = 14)

effect_indel_fig = plot_effect(indel_annotate)
ggsave(file.path(plotdir, "effect_indel.pdf"), effect_indel_fig)

# Plot mutated genes for snvs
# First determine length of genes
listMarts()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_info = getBM(attributes = c("start_position", "end_position", "strand", "chromosome_name", "hgnc_symbol"), filters = "chromosome_name", value = "MT", mart = ensembl)
somatic_variant_annotate_length = gene_info %>% 
  dplyr::mutate(length = end_position - start_position) %>% 
  dplyr::select(gene = hgnc_symbol, length, strand) %>% 
  right_join(somatic_variant_annotate, by = c("gene"))
saveRDS(somatic_variant_annotate_length, file.path(r_wd, "somatic_variant_annotate_length.rds"))

gene_fig = plot_gene(somatic_variant_annotate_length)
ggsave(file.path(plotdir, "gene.pdf"), gene_fig)

gene_length_fig = plot_gene(somatic_variant_annotate_length, length_correction = TRUE)
ggsave(file.path(plotdir, "gene_length_correct.pdf"), gene_length_fig)


### Modeling nr mutations

# Healthy age line
healthy_freq <- subset(total_sample_mutation_freq, 
                       (state == 'healthy') &
                         !(patient == 'AHH1') &
                         bulk == 'clone') %>% 
  dplyr::mutate(cov = factor(cov, levels = c(15, 30)))

glm_mixed_log_m <- glmer(freq ~ age + (0 + age | patient), 
                            data = healthy_freq, 
                            family = poisson(link = "log"))

glm_m = glm(freq ~ age, data = healthy_freq, family = poisson(link = "identity"))
glm_log_m = glm(freq ~ age, data = healthy_freq, family = poisson(link = "log"))


#Zero inflated model
library(pscl)
zeroinfl_m <- zeroinfl(freq ~ age, data = healthy_freq, dist="poisson")

# Alternatively fit random intercept
glm_mixed_log_randomintercept_m = glmer(freq ~ age + (1 | patient), 
      data = healthy_freq, 
      family = poisson(link = "log"))

# Linear mixed model
linear_mixed_m = lme(freq ~ age, random = ~-1 + age | patient, data = healthy_freq)
summary(linear_mixed_m)

# Simple linear model
linear_m = lm(freq ~ age, data = healthy_freq)

AIC(zeroinfl_m, glm_mixed_log_m, glm_mixed_log_randomintercept_m, glm_log_m, glm_m, linear_mixed_m, linear_m)
BIC(zeroinfl_m, glm_mixed_log_m, glm_mixed_log_randomintercept_m, glm_log_m, glm_m, linear_mixed_m, linear_m)

blood_glm_m = glm_m
get_ci_params(glm_m)
saveRDS(blood_glm_m, file.path(model_wd, "ageline.rds"))
blood_ageline_fig = plot_ageline(blood_glm_m)
ggsave(file.path(plotdir, "mut_accumulation", "blood_ageline.pdf"), blood_ageline_fig)
saveRDS(blood_ageline_fig, file.path(plotdir, "mut_accumulation", "blood_ageline.rds"))

genome_length = 2745186691
mt_genome_length = 16569
mt_rate_per_bp = coefficients(blood_glm_m)[2] / mt_genome_length
autosome_rate_per_bp = 14.2 / genome_length
mt_rate_per_bp / autosome_rate_per_bp

# Create colon model
colon_freq <- subset(total_sample_mutation_freq, 
                       (state == 'healthy_colon') &
                         !(patient == 'AHH1') &
                         bulk == 'clone')

colon_glm_m = glm(freq ~ age, data = colon_freq, family = poisson(link = "identity"))
saveRDS(colon_glm_m, file.path(model_wd, "colon_ageline.rds"))
colon_ageline_fig = plot_ageline(colon_glm_m)
ggsave(file.path(plotdir, "mut_accumulation", "colon_ageline.pdf"), colon_ageline_fig)
saveRDS(colon_ageline_fig, file.path(plotdir, "mut_accumulation", "colon_ageline.rds"))
get_ci_params(colon_glm_m)

# Create intestine model
intestine_freq <- subset(total_sample_mutation_freq, 
                     (state == 'healthy_intestine') &
                       !(patient == 'AHH1') &
                       bulk == 'clone')

intestine_glm_m = glm(freq ~ age, data = intestine_freq, family = poisson(link = "identity"))
saveRDS(intestine_glm_m, file.path(model_wd, "intestine_ageline.rds"))
intestine_ageline_fig = plot_ageline(intestine_glm_m)
ggsave(file.path(plotdir, "mut_accumulation", "intestine_ageline.pdf"), intestine_ageline_fig)
saveRDS(intestine_ageline_fig, file.path(plotdir, "mut_accumulation", "intestine_ageline.rds"))
get_ci_params(intestine_glm_m)

# Look at liver. Don't create a model, because of a lack of samples.
liver_freq <- subset(total_sample_mutation_freq, 
                     (state == 'healthy_liver') &
                       !(patient == 'AHH1') &
                       bulk == 'clone')
saveRDS(liver_freq, file.path(r_wd, "liver_freq.rds"))
liver_ageline_fig = plot_liver_ageline(liver_freq)
ggsave(file.path(plotdir, "mut_accumulation", "liver_ageline.pdf"), liver_ageline_fig)

# Tissues combined
tissue_freq = rbind(healthy_freq, colon_freq, intestine_freq) %>% 
  dplyr::mutate(state = factor(state, levels = c("healthy", "healthy_colon", "healthy_intestine")))
tissue_glm_m = glm(freq ~ age + state, data = tissue_freq, family = poisson(link = "identity"))
summary(tissue_glm_m)

# Does sequencing depth have an effect by itself? No.
cov_m = glm(freq ~ age + cov, data = healthy_freq, family = poisson(link = "identity"))
summary(cov_m)
saveRDS(cov_m, file.path(model_wd, "sequencing_depth.rds"))
seqdepth_fig = plot_cov_model2(cov_m)
ggsave(file.path(plotdir, "mut_accumulation", "sequencing_depth.pdf"), seqdepth_fig)

# Does HSCT have more mutations
hsct_freq = subset(total_sample_mutation_freq,
                   !(patient == "AHH1") &
                     (state == "healthy" | state == "recipient") &
                     bulk == "clone") %>% 
  dplyr::mutate(state = dplyr::recode(state, "healthy" = "Healthy blood", "recipient" = "HSCT recipient"),
                state = factor(state, levels = c("Healthy blood", "HSCT recipient")))
glm_m_mixed = glmer(freq ~ state + age + (1 | patient), family = poisson(link = "log"), data = hsct_freq)
glm_m = glm(freq ~ state + age, family = poisson(link = "identity"), data = hsct_freq)
BIC(glm_m_mixed, glm_m)
AIC(glm_m_mixed, glm_m)

saveRDS(glm_m, file.path(model_wd, "hsct_ageline.rds"))
hsct_ageline_fig = plot_hsct_model(glm_m)
ggsave(file.path(plotdir, "mut_accumulation", "hsct_ageline.pdf"), hsct_ageline_fig)

# Plot mean of hsct, to compare recipient and donor of same patient/hsct
mean_hsct_fig = plot_mean_hsct_freq(hsct_freq)
ggsave(file.path(plotdir, "mut_accumulation", "hsct_meanfreq.pdf"), mean_hsct_fig)


# Combine DX1, FU and DX2.
combi_dx1_fu_dx2_freq = dplyr::filter(total_sample_mutation_freq,
                                      (state == "DX1" | state == "FU" | state == "DX2") &
                                        bulk == "clone" &
                                        !(patient == "AHH1"))
m = glm(freq ~ state + age, family = poisson(link = "identity"), data = combi_dx1_fu_dx2_freq)
m2 = glm(freq ~ age, family = poisson(link = "identity"), data = combi_dx1_fu_dx2_freq)
BIC(m, m2)
AIC(m, m2)
anova(m2, m, test = "Chisq")

combi_dx1_fu_dx2_freq = dplyr::filter(total_sample_mutation_freq,
                                  (state == "healthy" | state == "DX1" | state == "FU" | state == "DX2") &
                                    bulk == "clone" &
                                    !(patient == "AHH1")) %>% 
  dplyr::mutate(state = dplyr::recode(state, "healthy" = "Healthy blood", "DX1" = "HSPCs leukemia patients", "FU" = "HSPCs leukemia patients", "DX2" = "HSPCs leukemia patients"),
                state = factor(state, levels = c("Healthy blood", "HSPCs leukemia patients")))
m = glm(freq ~ state + age, family = poisson(link = "identity"), data = combi_dx1_fu_dx2_freq)
saveRDS(m, file.path(model_wd, "chemo_dx1_fu_dx2_ageline.rds"))
dx1_fu_dx2_ageline_fig = plot_chemo_model(m)
ggsave(file.path(plotdir, "mut_accumulation", "dx1_fu_dx2_ageline.pdf"), dx1_fu_dx2_ageline_fig)

outlier_pvals = resid(m, type = "pearson") %>%
  abs() %>%
  multiply_by(-1) %>%
  pnorm() %>%
  multiply_by(2)
outlier_pvals_tb = tibble("sample" = combi_dx1_fu_dx2_freq$sample, "p" = outlier_pvals, 
                          "freq" = combi_dx1_fu_dx2_freq$freq, "state" = combi_dx1_fu_dx2_freq$state) %>%
  dplyr::filter(p < 0.05)
write_tsv(outlier_pvals_tb, file.path(r_wd, "dx1_fu_dx2_outliers.txt"))


# Do the chemo 2nd cancer samples have more mutations?
chemo_2ndcancer_freq = dplyr::filter(total_sample_mutation_freq,
                                  (state == "healthy" | chemo %in% c("DX2_cancer")) &
                                    bulk == "clone" &
                                    !(patient == "AHH1")) %>% 
  dplyr::mutate(state = dplyr::recode(state, "AML" = "2nd Leukemia", "ALL" = "2nd Leukemia", "healthy" = "Healthy blood"),
                state = factor(state, levels = c("Healthy blood", "2nd Leukemia")))
m = glm(freq ~ state + age, family = poisson(link = "identity"), data = chemo_2ndcancer_freq)
saveRDS(m, file.path(model_wd, "chemo_2nd_cancer_ageline.rds"))
chemo_2ndcancer_ageline_fig = plot_chemo_model(m)
ggsave(file.path(plotdir, "mut_accumulation", "chemo_2nd_cancer_ageline.pdf"), chemo_2ndcancer_ageline_fig)


outlier_pvals = resid(m, type = "pearson") %>%
  abs() %>%
  multiply_by(-1) %>%
  pnorm() %>%
  multiply_by(2)
outlier_pvals_tb = tibble("sample" = chemo_2ndcancer_freq$sample, "p" = outlier_pvals, 
                          "freq" = chemo_2ndcancer_freq$freq, "state" = chemo_2ndcancer_freq$state) %>%
  dplyr::filter(p < 0.05)
write_tsv(outlier_pvals_tb, file.path(r_wd, "chemo_2nd_cancer_outliers.txt"))


#Do the CB samples treated with different chemos show differences.
cb_chemo_freq = dplyr::filter(total_sample_mutation_freq, state == "CB_Chemo") %>% 
  dplyr::mutate(chemo = stringr::str_remove(chemo, "CB_"),
                chemo = factor(chemo, levels = c("CTRL", 
                                                 "CAR", 
                                                 "CISPL", 
                                                 "CYTA", 
                                                 "DOX", 
                                                 "MAPH", 
                                                 "RAD", 
                                                 "VINCRIS", 
                                                 "GCV",
                                                 "FC",
                                                 "GCV+FC")))
saveRDS(cb_chemo_freq, file.path(r_wd, "CB_chemo.rds"))

cb_chemo_fig = plot_cb_chemo(cb_chemo_freq)
ggsave(file.path(plotdir, "mut_accumulation", "CB_chemo.pdf"), cb_chemo_fig)


### Modeling of CNV
# ageline
healthy_cnv <- as.data.frame(subset(cnv_df, 
                      (state == 'healthy') &
                        !(patient == 'AHH1') &
                        mt_mean > 1000 &
                        bulk == 'clone'))

ageline_cnv_model <- lme(cnv_mean ~ age, 
                         random = ~ 0 + age | patient, 
                         data = healthy_cnv)
ageline_cnv_model2 <- lme(cnv_mean ~ age, 
                         random = ~ 1 | patient, 
                         data = healthy_cnv)

cnv_lm_m = lm(cnv_mean ~ age, data = healthy_cnv)

BIC(ageline_cnv_model, ageline_cnv_model2, cnv_lm_m)
AIC(ageline_cnv_model, ageline_cnv_model2, cnv_lm_m)

saveRDS(ageline_cnv_model, file.path(model_wd, "cnv_ageline.rds"))
cnv_ageline_fig = plot_cnv_ageline(ageline_cnv_model)
ggsave(file.path(plotdir, "copy_numbers", "cnv_ageline.pdf"), cnv_ageline_fig)

colon_cnv <- as.data.frame(subset(cnv_df, 
                      (state == 'healthy_colon') &
                        !(patient == 'AHH1') &
                        mt_mean > 1000 &
                        bulk == 'clone'))

colon_cnv_model <- lme(cnv_mean ~ age, 
                         random = ~ 0 + age | patient, 
                         data = colon_cnv)

saveRDS(colon_cnv_model, file.path(model_wd, "colon_cnv_ageline.rds"))
cnv_ageline_fig = plot_cnv_ageline(colon_cnv_model)
ggsave(file.path(plotdir, "copy_numbers", "colon_cnv_ageline.pdf"), cnv_ageline_fig)

intestine_cnv <- as.data.frame(subset(cnv_df, 
                    (state == 'healthy_intestine') &
                      !(patient == 'AHH1') &
                      mt_mean > 1000 &
                      bulk == 'clone'))

intestine_cnv_model <- lme(cnv_mean ~ age, 
                       random = ~ 0 + age | patient, 
                       data = intestine_cnv)

saveRDS(intestine_cnv_model, file.path(model_wd, "intestine_cnv_ageline.rds"))
cnv_ageline_fig = plot_cnv_ageline(intestine_cnv_model)
ggsave(file.path(plotdir, "copy_numbers", "intestine_cnv_ageline.pdf"), cnv_ageline_fig)

liver_cnv <- as.data.frame(subset(cnv_df, 
                    (state == 'healthy_liver') &
                      !(patient == 'AHH1') &
                      mt_mean > 1000 &
                      bulk == 'clone'))
saveRDS(liver_cnv, file.path(r_wd, "liver_cnv.rds"))
cnv_ageline_fig = plot_cnv_liver_ageline(liver_cnv)
ggsave(file.path(plotdir, "copy_numbers", "liver_cnv_ageline.pdf"), cnv_ageline_fig)


# Compare mean cnv between tissues
tissue_cnv = rbind(healthy_cnv, colon_cnv, intestine_cnv, liver_cnv)
tissue_model = lme(cnv_mean ~ state, 
    random = ~ 1 | patient, 
    data = tissue_cnv)
summary(tissue_model)
ggpredict(tissue_model)

colon_intestine_cnv = rbind(colon_cnv, intestine_cnv)
colon_intestine_model = lme(cnv_mean ~ state, 
                   random = ~ 1 | patient, 
                   data = colon_intestine_cnv)
summary(colon_intestine_model)
ggpredict(colon_intestine_model)


# Check if age has an effect if all tissues are pooled (No)
combi_tissues_cnv_model <- lme(cnv_mean ~ age + state, 
                               random = ~ 0 + age | patient, 
                               data = tissue_cnv)


# Does cov have an effect by itself:
cov_m <- lme(cnv_mean ~ cov, random = ~ 1 | patient, 
               data = foetus_cnv)
summary(cov_m)

saveRDS(cov_m, file.path(model_wd, "cnv_fetus_cov.rds"))
cov_fig = plot_cnv_model(cov_m, var = "cov")
ggsave(file.path(plotdir, "copy_numbers", "cnv_fetus_cov.pdf"), cov_fig)


cov_m <- lme(cnv_mean ~ cov, random = ~ 1 | patient, 
             data = healthy_cnv)
summary(cov_m)

saveRDS(cov_m, file.path(model_wd, "cnv_cov.rds"))
cov_fig = plot_cnv_model(cov_m, var = "cov")
ggsave(file.path(plotdir, "copy_numbers", "cnv_cov.pdf"), cov_fig)


# Bulk comparison. No difference. Model with different variance between groups also not significant.
clone_bulk_cnv <- as.data.frame(subset(cnv_df, 
                          state == 'healthy' &
                            !(patient == 'AHH1') &
                           mt_mean > 1000 &
                            source == "Boxtel")) %>% 
  dplyr::mutate(bulk = factor(bulk, levels = c("clone", "bulk")))
bulk_m = lme(cnv_mean ~ age + bulk, random = ~ 1 | patient, data = clone_bulk_cnv)
bulk_m2 = lme(cnv_mean ~ bulk, random = ~ 1 | patient, data = clone_bulk_cnv)
bulk_m3 = lme(cnv_mean ~ bulk, random = ~ 1 | patient, data = clone_bulk_cnv, weights = varIdent(form = ~1|bulk))
BIC(bulk_m, bulk_m2, bulk_m3)
AIC(bulk_m, bulk_m2, bulk_m3)
summary(bulk_m2)
saveRDS(bulk_m2, file.path(model_wd, "cnv_bulk.rds"))
bulk_fig = plot_cnv_model(bulk_m2, var = "bulk")
ggsave(file.path(plotdir, "copy_numbers", "cnv_bulk.pdf"), bulk_fig)

# AML is significant. Highly correlated with source. Source is even more significant.
aml_cnv = as.data.frame(subset(cnv_df,
                  !(patient == "AHH1") &
                   (state == "AML" & chemo == "DX1_cancer") &
                   mt_mean > 1000 &
                    bulk == "clone")) %>% 
  dplyr::mutate(state = factor(state, levels = unique(state)),
                source = dplyr::recode(source, "Boxtel" = "In-house"))



saveRDS(aml_cnv, file.path(r_wd, "aml_cnv.rds"))
aml_cnv_fig = plot_aml_cnv_model(aml_cnv)
ggsave(file.path(plotdir, "copy_numbers", "cnv_aml.pdf"), aml_cnv_fig)

# Does hsct have an effect on cnv
hsct_cnv = subset(cnv_df,
                  !(patient == "AHH1") &
                    (state == "healthy" | state == "recipient") &
                    mt_mean > 1000 &
                    bulk == "clone" &
                    source == "Boxtel") %>% 
  dplyr::mutate(state = dplyr::recode(state, "healthy" = "Healthy blood", "recipient" = "HSCT recipient"),
                state = factor(state, levels = c("Healthy blood", "HSCT recipient")))


hsct_m = lme(cnv_mean ~ age + state, random = ~ 1 | patient, data = hsct_cnv)
hsct_m2 = lme(cnv_mean ~ state, random = ~ 1 | patient, data = hsct_cnv)
BIC(hsct_m, hsct_m2)
AIC(hsct_m, hsct_m2)



saveRDS(hsct_m2, file.path(model_wd, "cnv_hsct.rds"))
hsct_cnv_fig = plot_cnv_model(hsct_m2)
ggsave(file.path(plotdir, "copy_numbers", "cnv_hsct.pdf"), hsct_cnv_fig)

# Plot mean of hsct, to compare recipient and donor of same patient/hsct
mean_hsct_cnv_fig = plot_mean_hsct_cnv(hsct_cnv)
ggsave(file.path(plotdir, "copy_numbers", "hsct_meancnv.pdf"), mean_hsct_cnv_fig)


# Check if there is an effect of the time after transplant
time_after_transplant = tibble("patient" = c("HSCT2", "HSCT3", "HSCT4", "HSCT8", "HSCT10", "9386", "HSCT13", "HSCT14", "HSCT5"), 
                               "time_transplant" = c(0.416666667, 1.416666667, 3.333333333, 0.083333333, 2.416666667, 0.25, 24.5833, 16.25, 2.5))
recipient_cnv = hsct_cnv %>% 
  dplyr::filter(hsct == "recipient") %>% 
  left_join(time_after_transplant, by = "patient")
recipient_m = lme(cnv_mean ~ time_transplant, random = ~ 1 | patient, data = recipient_cnv)
saveRDS(recipient_m, file.path(model_wd, "cnv_hsct_time.rds"))
recipient_fig = plot_recipient_time_model(recipient_m)
ggsave(file.path(plotdir, "copy_numbers", "hsct_time_transplant.pdf"), recipient_fig)

hsct_sig_samples = c("9386FU2D19", "9386FUHSC1B15", "9386FUHSC1P23", "PMC9386-DX2BMWT-HSP1C8", "PMC9386-DX2BMWT-HSP2I16", 
                     "PMC9386-DX2BMWT-HSP2I8", "FPHSCT10RPBMPP2", "FPHSCT10RPBMPP8", "PMCHSCT10-RPBWT-HSPHSCT10-2J12", "PMCHSCT10-RPBWT-HSPHSCT10-1J15")

recipient_cnv = recipient_cnv %>% 
  dplyr::mutate(hsct_sig = ifelse(sample %in% hsct_sig_samples, "SBSA", "not_SBSA"))

# Check if the blood type influences copy numbers.
sample_type_eline = read_tsv("samples_types_Eline.txt", col_types = c("ccc"))
pb_bm_cnv_df = subset(cnv_df,
                  !(patient == "AHH1") &
                    (state == "healthy" | state == "recipient" | state == "DX1" | state == "DX2" | state == "FU") &
                    mt_mean > 1000 &
                    bulk == "clone" &
                    source == "Boxtel") %>% 
  dplyr::left_join(sample_type_eline, by = c("patient", "state")) %>% 
  dplyr::mutate(state = dplyr::recode(state, "healthy" = "Normal blood", "recipient" = "HSCT recipient", "DX1" = "Normal blood", 
                                      "DX2" = "Normal blood", "FU" = "Normal blood"),
                state = factor(state, levels = c("Normal blood", "HSCT recipient"))) %>% 
  dplyr::mutate(blood_type = case_when(
    sample == "N01BMHSPCCB8" ~ "BM",
    patient == "CB112" ~ "CB",
    patient %in% c("HSCT13", "HSCT14") ~ "PB",
    patient %in% c("HSCT8") ~ "BM",
    patient %in% c("AC33", "AC41", "ACC55", "AC63", "BCH") ~ "BM",
    patient %in% c("PMC16332", "PMC17556", "PMC20106", "PMC21586", "PMC21636", "PMC22813", "PMC07276", "PMC09357") ~ "BM",
    patient %in% c("MH2", "NR1", "NR2") ~ "Li",
    hsct == "donor" ~ "BM",
    hsct == "recipient" ~ "PB"
  )) %>% 
  dplyr::mutate(blood_type = ifelse(is.na(Sample_type), blood_type, Sample_type)) %>% 
  dplyr::filter(!is.na(blood_type))


blood_type_names = tibble("blood_type" = c("BM", "CB", "Li", "PB"), "full_blood_type" = c("Bone marrow", "Cord blood", "Liver", "Peripheral blood"))
pb_bm_cnv_df = dplyr::left_join(pb_bm_cnv_df, blood_type_names, by = "blood_type")
bloodtype_m = lme(cnv_mean ~ state + blood_type, random = ~ 1 | patient, data = pb_bm_cnv_df)
saveRDS(bloodtype_m, file.path(model_wd, "cnv_bloodtype_transplant.rds"))
bloodtype_cnv_fig = plot_cnv_model(bloodtype_m, var = "full_blood_type", col = "state", remove_guide = FALSE, plot_p = FALSE)
bloodtype_cnv_fig2 = plot_cnv_model(bloodtype_m, var = "state", col = "blood_type", remove_guide = FALSE, plot_p = FALSE)
ggsave(file.path(plotdir, "copy_numbers", "hsct_bloodtype2.pdf"), ggarrange(bloodtype_cnv_fig, bloodtype_cnv_fig2))
ggsave(file.path(plotdir, "copy_numbers", "hsct_bloodtype3.pdf"), bloodtype_cnv_fig)

# effect of blood type
bloodtype_m2 = lme(cnv_mean ~ blood_type, random = ~ 1 | patient, data = pb_bm_cnv_df)
saveRDS(bloodtype_m2, file.path(model_wd, "cnv_bloodtype.rds"))

# Restrict to non-transplant
nontrans_cnv_df = dplyr::filter(pb_bm_cnv_df, state == "Normal blood")
summary(lme(cnv_mean ~ blood_type, random = ~ 1 | patient, data = nontrans_cnv_df))

# Effect of transplant or not
bloodtype_m3 = lme(cnv_mean ~ state, random = ~ 1 | patient, data = pb_bm_cnv_df)
saveRDS(bloodtype_m3, file.path(model_wd, "cnv_transplant.rds"))

# Restrict to PB
pb_cnv_df = dplyr::filter(pb_bm_cnv_df, blood_type == "PB")
summary(lme(cnv_mean ~ state, random = ~ 1 | patient, data = pb_cnv_df))

# Make figure with only BM and PB
pb_bm_cnv_df2 = dplyr::filter(pb_bm_cnv_df, blood_type %in% c("BM", "PB"))
bloodtype_m = lme(cnv_mean ~ blood_type, random = ~ 1 | patient, data = pb_bm_cnv_df2)
bloodtype_cnv_fig = plot_cnv_model(bloodtype_m, var = "blood_type", col = "state", remove_guide = FALSE, plot_p = T)
ggsave(file.path(plotdir, "copy_numbers", "hsct_bloodtype4.pdf"), bloodtype_cnv_fig)


# Compare all the cnvs of the chemo/2nd cancer project
chemo_cnv = dplyr::mutate(cnv_df, state = ifelse(chemo == "DX1_cancer" | chemo == "DX2_cancer", chemo, state)) %>% 
  dplyr::filter(state %in% c("healthy", "DX1", "FU", "DX2", "DX1_cancer", "DX2_cancer") &
                  bulk == "clone" &
                  source == "Boxtel" &
                  !(patient == "AHH1")) %>% 
  dplyr::mutate(state = dplyr::recode(state, "healthy" = "Healthy blood", 
                                 "DX1" = "Diagnosis", 
                                 "DX1_cancer" = "Leukemia", 
                                 "FU" = "Follow-up", 
                                 "DX2" = "Diagnosis 2", 
                                 "DX2_cancer" = "2nd Leukemia"),
                state = factor(state, levels = c("Healthy blood", "Diagnosis", "Leukemia", "Follow-up", "Diagnosis 2", "2nd Leukemia")))
m = lme(cnv_mean ~ state + age, random =  ~1 | patient, data = chemo_cnv)
m2 = lme(cnv_mean ~ state, random =  ~1 | patient, data = chemo_cnv)
BIC(m, m2)
AIC(m, m2)
saveRDS(m2, file.path(model_wd, "cnv_chemo_2ndcancer.rds"))
chemo_cnv_fig = plot_cnv_model(m2, plot_p = FALSE)
ggsave(file.path(plotdir, "copy_numbers", "chemo_2ndcancer_cnv.pdf"), chemo_cnv_fig)

sum_m <- summary(m2)
pvals <- formatC(sum_m$tTable[2:6,5],
                digits = 3,
                format = 'f')
padj = p.adjust(pvals, method = "fdr")
p_tb = tibble("name" = names(pvals), "pval" = pvals, "padj" = padj)
write_tsv(p_tb, file.path(r_wd, "cnv_chemo_2ndcancer.txt"))

# Find the outliers
outlier_pvals = resid(m2, type = "pearson") %>% 
  abs() %>%
  multiply_by(-1) %>%
  pnorm() %>%
  multiply_by(2)
outlier_pvals_tb = tibble("sample" = chemo_cnv$sample, "p" = outlier_pvals, "state" = chemo_cnv$state, "cnv" = chemo_cnv$cnv_mean,) %>% 
  dplyr::filter(p < 0.05)
write_tsv(outlier_pvals_tb, file.path(r_wd, "cnv_chemo_2ndcancer_outliers.txt"))
outlier_plot = plot(m2, id = 0.05, idLabels = chemo_cnv$sample)
pdf(file.path(plotdir, "copy_numbers", "chemo_2ndcancer_cnv_outliers.pdf"))
print(outlier_plot)
dev.off()

variance_tb = broom::tidy(fligner.test(cnv_mean ~ state, data = chemo_cnv))
write_tsv(variance_tb, file.path(r_wd, "cnv_chemo_variance_test.txt"))

#Do the CB samples treated with different chemos show differences.
cb_chemo_cnv = dplyr::filter(cnv_df, state == "CB_Chemo") %>% 
  dplyr::mutate(chemo = stringr::str_remove(chemo, "CB_"),
                chemo = factor(chemo, levels = c("CTRL", 
                                                 "CAR", 
                                                 "CISPL", 
                                                 "CYTA", 
                                                 "DOX", 
                                                 "MAPH", 
                                                 "RAD", 
                                                 "VINCRIS", 
                                                 "GCV",
                                                 "FC",
                                                 "GCV+FC")))
saveRDS(cb_chemo_cnv, file.path(r_wd, "CB_chemo_cnv.rds"))

cb_chemo_fig = plot_cb_chemo_cnv(cb_chemo_cnv)
ggsave(file.path(plotdir, "copy_numbers", "CB_chemo_cnv.pdf"), cb_chemo_fig)



#Compare copy number between bulk and clone in TARGET data
target_cnv = cnv_df %>% dplyr::filter(source == "TARGET")
saveRDS(target_cnv, file.path(r_wd, "target_cnv.rds"))
aml_bulk_cnv_fig = compare_aml_bulk_cnv(target_cnv)
ggsave(file.path(plotdir, "copy_numbers", "aml_bulk_cnv_comp.pdf"), aml_bulk_cnv_fig)

# This difference is even stronger when looking at only mt reads.
# It is thus likely not caused by a decreased ploidy of the aml.
aml_bulk_mtreads_fig = compare_aml_bulk_mt_reads(target_cnv)
ggsave(file.path(plotdir, "copy_numbers", "aml_bulk_mt_reads.pdf"), aml_bulk_mtreads_fig)
broom::tidy(wilcox.test(n_mean ~ bulk, data = target_cnv, paired = TRUE))

# Check if there is a relation between the number of mutations and the copy numbers
healthy_combi = left_join(healthy_freq, healthy_cnv)
m = glm(freq ~ cnv_mean, data = healthy_combi, family = poisson(link = "identity"))
m2 = lm(freq ~ cnv_mean, data = healthy_combi)
AIC(m, m2)
BIC(m, m2)
saveRDS(m, file.path(model_wd, "freq_cnv.rds"))
summary(m)$coefficients
freq_cnv_fig = plot_freq_cnv_model(m)
ggsave(file.path(plotdir, "healthy_freq_cnv_correlation.pdf"), freq_cnv_fig)


# Check if the samples with more nuclear mutations also have more mitochondrial mutations

# Load nuclear data
fnames = list.files(path="~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/AgeLine/hg38/SNV/",pattern="ageline_table.txt",recursive = T, full.names = T)
nuclear_muts = purrr::map(fnames, read_delim, delim = " ", col_type = c("cddcdccdl")) %>% 
  bind_rows() %>% 
  dplyr::select(sample = sample_id, norm_muts)

#Determine samples that are above the line in mito
healthy_freq$high_mut = predict(blood_glm_m) > healthy_freq$freq

# Combine nuclear and mito data
nuclear_mito_muts = inner_join(nuclear_muts, healthy_freq, by = "sample")
cor(nuclear_mito_muts$norm_muts, nuclear_mito_muts$freq)
ggplot(nuclear_mito_muts, aes(x = freq, y = norm_muts)) +
  geom_point()

# freq is not useful in model
m = lme(norm_muts ~ age, random = ~ -1 + age | patient, data = nuclear_mito_muts)
# m2 = lme(norm_muts ~ age + freq, random = ~ -1 + age | patient, data = nuclear_mito_muts)
# summary(m2)
# vif(m2)

# Check which samples are above the ageline in nuclear
nuclear_mito_muts$nuclear_high_mut = predict(m) > nuclear_mito_muts$norm_muts

# Check if the samples that are above the ageline in nuclear are the same as the samples in mito
higher_lower_muts = nuclear_mito_muts %>% 
  dplyr::select(high_mut, nuclear_high_mut) %>% 
  dplyr::mutate(high_mut = ifelse(high_mut, "Above", "Below"),
                high_mut = factor(high_mut, levels = c("Below", "Above")),
                nuclear_high_mut = ifelse(nuclear_high_mut, "Above", "Below"),
                nuclear_high_mut = factor(nuclear_high_mut, levels = c("Below", "Above"))) %>% 
  dplyr::group_by(high_mut, nuclear_high_mut) %>% 
  dplyr::count() %>% 
  dplyr::ungroup()
saveRDS(higher_lower_muts, file.path(r_wd, "higher_lower_muts.rds"))

chisq_res = higher_lower_muts %>% 
  tidyr::pivot_wider(values_from = n, names_from = nuclear_high_mut) %>% 
  dplyr::select(-high_mut) %>% 
  chisq.test(simulate.p.value = TRUE) %>% 
  broom::tidy()

above_below_ageline_fig = plot_above_below_ageline(higher_lower_muts)
ggsave(file.path(plotdir, "above_below_agelines_nuclear_mito.pdf"), above_below_ageline_fig)


# Look at pcawg data.
pcawg_cancer_type = read_tsv("PCAWG_donors_metadata.tsv", col_types = cols_only(submitted_donor_id = col_character(),
                                                                                histology_abbreviation = col_character(),
                                                                                donor_age_at_diagnosis = col_character()))

# Read the pcawg file with the number of mito mutations per sample
mito_count = read_tsv("samples_with_mtDNA_mutations_forFreekManders.txt") %>% 
  dplyr::mutate(temp = submitter_donor_id) %>% 
  tidyr::separate(temp, into = c("something", "sample_id"), sep = "::") %>% 
  dplyr::select(-something)

# Combine the mutation counts with meta data
pcawg = inner_join(mito_count, pcawg_cancer_type, by = c("sample_id" = "submitted_donor_id"))


# PCAWG frequency amounts
pcawg_freq = pcawg %>%
  dplyr::select(type = histology_abbreviation, freq = nr_snvs)

blood_cancer_types = c("Lymph-BNHL", 
                       "Lymph-CLL", 
                       "Lymph-NOS", 
                       "Myeloid-AML", 
                       "Myeloid-AML,Myeloid-MPN", 
                       "Myeloid-MDS", 
                       "Myeloid-MPN",
                       "ALL",
                       "AML")
colon_cancer_types = c("ColoRect-AdenoCA")
####intestine_cancer_types = c("ColoRect−AdenoCA") CHECK THIS. IS THERE SOMETHING I CAN USE
liver_cancer_types = c("Biliary-AdenoCA")
blood_high_mut_samples = compare_pcawg_freq(pcawg_freq, "healthy", "Healthy blood", total_sample_mutation_freq, blood_cancer_types, "Blood cancer", blood_glm_m, pcawg)
colon_high_mut_samples = compare_pcawg_freq(pcawg_freq, "healthy_colon", "Normal colon", total_sample_mutation_freq, colon_cancer_types, "Colon cancer", colon_glm_m, pcawg)
intestine_high_mut_samples = compare_pcawg_freq(pcawg_freq, "healthy_intestine", "Normal intestine", total_sample_mutation_freq, intestine_cancer_types, "small_intestine_cancer", intestine_glm_m, pcawg)
liver_high_mut_samples = compare_pcawg_freq(pcawg_freq, "healthy_liver", "Normal liver", total_sample_mutation_freq, liver_cancer_types, "liver_cancer", liver_glm_m, pcawg)

pcawg_freq_blood = dplyr::filter(pcawg_f)


# pcawg spectra
pcawg_snv = read_tsv("TCMA-MutationSNV.tsv") %>% 
  dplyr::rename(mut_pos = position, alt = var)

cancer_snv = cancer_variants %>%
  dplyr::mutate(chrom = "MT") %>% 
  dplyr::select(sample_id = sample, cancer_type = state, chrom, mut_pos, ref, alt, var_type = effect_type)

pcawg_snv_incl_our = rbind(pcawg_snv, cancer_snv)


plot_save_pcawg_spectra_tissue(pcawg_snv_incl_our, healthy_variants, "Healthy blood", "Blood cancer", blood_high_mut_samples, blood_cancer_types)
plot_save_pcawg_spectra_tissue(pcawg_snv, colon_variants, "Normal colon", "Colon cancer", colon_high_mut_samples, colon_cancer_types)
plot_save_pcawg_spectra_tissue(pcawg_snv, intestine_variants, "Normal intestine", "pcawg_intestine", intestine_high_mut_samples, intestine_cancer_types)
plot_save_pcawg_spectra_tissue(pcawg_snv, liver_variants, "Normal liver", "Liver cancer", liver_high_mut_samples, liver_cancer_types)

# Compare light strand of blood cancer with colon cancer.
blood_pcawg_snv = pcawg_snv_incl_our %>% 
  dplyr::filter(cancer_type %in% blood_cancer_types)
colon_pcawg_snv = pcawg_snv %>% 
  dplyr::filter(cancer_type %in% colon_cancer_types)
spec_blood_l = get_mt_split_spectra(blood_pcawg_snv, ref_genome, remove_duplis = FALSE)
spec_colon_l = get_mt_split_spectra(colon_pcawg_snv, ref_genome, remove_duplis= FALSE)
spec_m = rbind(spec_blood_l$spec_c_t, spec_colon_l$spec_c_t)[,-3]
chisq.test(spec_m, simulate.p.value = TRUE)
cos_sim(unlist(spec_blood_l$spec_c_t), unlist(spec_colon_l$spec_c_t))


blood_gr = to_granges_mt(blood_pcawg_snv, remove_duplis = FALSE)
blood_gr = blood_gr[blood_gr$ref %in% c("C", "T")]
colon_gr = to_granges_mt(colon_pcawg_snv, remove_duplis = FALSE)
colon_gr = colon_gr[colon_gr$ref %in% c("C", "T")]

grl = GRangesList(blood_gr, colon_gr)
names(grl) = c(paste0("Blood cancer \nNo. mutations = ", length(blood_gr)),
               paste0("Colon cancer \nNo. mutations = ", length(colon_gr)))
mut_mat = mut_matrix(grl, ref_genome)
pcawg_profile_comparison_fig = plot_compare_profiles(mut_mat[,1], mut_mat[,2], profile_names = names(grl))
blood_colon_profile_fig = plot_96_profile(mut_mat)
ggsave(file.path(plotdir, "pcawg", "blood_colon_profile.pdf"), blood_colon_profile_fig)
saveRDS(mut_mat, file.path(plotdir, "pcawg", "blood_colon_mut_mat.rds"))



cos_sim(mut_mat[33:48,1], mut_mat[33:48,2])




# Determine if high mut pcawg samples have different nuclear sigs.
# blood_high_mut_samples = pcawg %>% dplyr::filter(nr_snvs >= 5 & histology_abbreviation %in% blood_cancer_types) %>% pull(sample_id)
# colon_high_mut_samples = pcawg %>% dplyr::filter(nr_snvs >= 5 & histology_abbreviation %in% colon_cancer_types) %>% pull(sample_id)

pcawg_meta_names = read_tsv("PCAWG_donors_metadata_exploded_ann_rep.tsv") %>% 
  dplyr::select(submitted_donor_id = "#submitted_donor_id", tumor_wgs_icgc_specimen_id, tumor_wgs_aliquot_id, tumor_rna_seq_aliquot_id)

pcawg_nuclear_sigs = read_csv("PCAWG_sigProfiler_SBS_signatures_in_samples.csv") %>% 
  left_join(pcawg_meta_names, by = c("Sample Names" = "tumor_wgs_icgc_specimen_id")) %>% 
  dplyr::select(-`Sample Names`, -tumor_wgs_aliquot_id, -tumor_rna_seq_aliquot_id)

create_pcawg_nuclear_sigs(pcawg_nuclear_sigs, blood_high_mut_samples, "blood")
create_pcawg_nuclear_sigs(pcawg_nuclear_sigs, colon_high_mut_samples, "colon")

# Determine if high mut pcawg samples have different drivers
# Read pcawg drivers and get correct names.
pcawg_drivers = read_tsv("TableS3_panorama_driver_mutations_ICGC_samples.public.tsv") %>% 
  left_join(pcawg_meta_names, by = c("sample_id" = "tumor_wgs_aliquot_id")) %>% 
  dplyr::select(-sample_id, -tumor_wgs_icgc_specimen_id, -tumor_rna_seq_aliquot_id) %>% 
  dplyr::filter(submitted_donor_id %in% pcawg$sample_id)

create_pcawg_drivers(pcawg_drivers, blood_cancer_types, blood_high_mut_samples, "blood")

# Determine if high mut pcawg samples have different expression.
expression_tbl = read_tsv("pcawg.rnaseq.transcript.expr.counts.tsv", col_names = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("Feature")

# Transform names and select only pcawg samples
samples_rna = colnames(expression_tbl)
samples_rna = tibble("sample_names" = samples_rna) %>% 
  left_join(pcawg_meta_names, by = c("sample_names" = "tumor_rna_seq_aliquot_id")) %>% 
  dplyr::select(-tumor_wgs_aliquot_id, -tumor_wgs_icgc_specimen_id, -sample_names) %>% 
  dplyr::pull(submitted_donor_id)
colnames(expression_tbl) = samples_rna
expression_tbl = expression_tbl[, samples_rna %in% pcawg$sample_id]

do_pcawg_expression_analysis(pcawg, expression_tbl, blood_cancer_types, blood_high_mut_samples, "blood")
do_pcawg_expression_analysis(pcawg, expression_tbl, colon_cancer_types, colon_high_mut_samples, "colon")

do_pcawg_expression_analysis = function(pcawg, expression_tbl, cancer_types, high_mut_samples, name){

# Select only tissue specific samples
tissue_samples = pcawg %>% 
  dplyr::filter(histology_abbreviation %in% cancer_types) %>% 
  pull(sample_id)
expression_tbl = expression_tbl[, colnames(expression_tbl) %in% tissue_samples]

# Create a matrix for DESEQ2
expression_m = expression_tbl %>% 
  round() %>% 
  as.matrix()
mode(expression_m) = "integer"

# Create metadata for DESEQ2
coldata = data.frame("condition" = factor(ifelse(colnames(expression_m) %in% high_mut_samples, "High_mut", "Other")))
rownames(coldata) = colnames(expression_m)

#Create DESEQ2 object and run analysis
dds <- DESeqDataSetFromMatrix(countData = expression_m,
                              colData = coldata,
                              design = ~ condition)


keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("Other","High_mut"))

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)


# Plot and write results
pdf(file.path(plotdir, "pcawg", paste0(name, "_diffexpr_unshrunk.pdf")))
plotMA(res, ylim=c(-3,3))
dev.off()

summary(res)
resOrdered <- res[order(res$pvalue),]
write_tsv(as.data.frame(resOrdered), file.path(r_wd, "pcawg_diffexpr", paste0(name, "_res.txt")))

# Go from transcripts to gene symbols
sig_transcripts = resOrdered %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < 0.05) %>% rownames()
expression_m[rownames(expression_m) %in% sig_transcripts ,coldata$condition == "High_mut"]
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
bm = getBM(c("ensembl_transcript_id", "ensembl_transcript_id_version", "hgnc_symbol"), mart = mart)
sig_transcripts_tbl = tibble("ensembl_transcript_id" = sig_transcripts) %>% 
  dplyr::mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.[0-9]"))
gene_symbols = left_join(sig_transcripts_tbl, bm, by = "ensembl_transcript_id") %>% 
  dplyr::select(hgnc_symbol) %>% 
  dplyr::filter(!duplicated(hgnc_symbol) & !is.na(hgnc_symbol))
write_tsv(gene_symbols, file.path(r_wd, "pcawg_diffexpr", paste0(name, "_genesymbols.txt")))

gene_symbols_all = tibble("ensembl_transcript_id" = rownames(dds)) %>% 
  dplyr::mutate(ensembl_transcript_id = str_remove(ensembl_transcript_id, "\\.[0-9]")) %>%
  left_join(bm, by = "ensembl_transcript_id") %>% 
  dplyr::select(hgnc_symbol) %>% 
  dplyr::filter(!duplicated(hgnc_symbol) & !is.na(hgnc_symbol))
write_tsv(gene_symbols_all, file.path(r_wd, "pcawg_diffexpr", paste0(name, "_genesymbols_all.txt")))

# Shrink fold changes
resLFC <- lfcShrink(dds, coef="condition_High_mut_vs_Other", type="apeglm")
write_tsv(as.data.frame(resLFC), file.path(r_wd, "pcawg_diffexpr", paste0(name, "_resLFC.txt")))
pdf(file.path(plotdir, "pcawg", paste0(name, "_diffexpr_shrunk.pdf")))
plotMA(resLFC, colLine = NA, ylim=c(-3,3))
dev.off()

}

# Compare pcawg copy number to healthy
pcawg_copy_number = read_tsv("TCMA-CopyNumber.tsv") %>% 
  dplyr::select(type = cancer_type, copy_number = tumor_copy_number, patient = sample_id)

compare_pcawg_cnv(pcawg_copy_number, cnv_df, "healthy", "Healthy blood", "Blood cancer", blood_cancer_types)
compare_pcawg_cnv(pcawg_copy_number, cnv_df, "healthy_colon", "Normal colon", "Colon cancer", colon_cancer_types)
compare_pcawg_cnv(pcawg_copy_number, cnv_df, "healthy_intestine", "Normal intestine", "Intestine cancer", intestine_cancer_types)
compare_pcawg_cnv(pcawg_copy_number, cnv_df, "healthy_liver", "Normal liver", "Liver scancer", liver_cancer_types)
