# Load VCFs and filter variants that PASSED previous filtering

read_all_vcfs <- function(wd, vcf_filter, sample_reference_df) {
  
  old_dir = setwd(wd)
  on.exit(setwd(old_dir), add = TRUE)
  
  # Load VCFs to analyse
  # Merged VCFs of same sample from Mito pipeline WDL
  
  # Find vcf files
  vcf_fnames <- list.files(pattern = vcf_filter, 
                           recursive = T)
  
  # Remove samples not in sample reference
  sample_names = vcf_fnames %>% 
    basename() %>% 
    str_remove("\\..*")
  sample_f = sample_names %in% sample_reference_df$sample
  vcf_fnames = vcf_fnames[sample_f]
  sample_names = sample_names[sample_f]
  
  # Check that no vcfs are absent
  samples_no_vcf = sample_reference_df$sample[!sample_reference_df$sample %in% sample_names]
  if (length(samples_no_vcf)){
    stop(paste0("These samples don't have a vcf: ", paste0(samples_no_vcf, collapse = ", ")))
  }
  
  
  # Not sure if this can't be removed
  vcf_fnames <- grep(pattern = "merged", 
                     x = vcf_fnames, 
                     val = T, 
                     invert = T)
  
  # Read data
  vcf_l <- lapply(vcf_fnames, readVcf)
  
  # Add sample names
  names(vcf_l) <- sample_names
  return(vcf_l)
}
filter_passed_variants <- function(vcf_l) {
  
  vcf_l = purrr::map(vcf_l, ~.x[granges(.x)$FILTER == "PASS"])
  
  return(vcf_l)
}

# Return all called variants in a DF

grab_all_variants<- function(vcf_l, mito_reference_df) {
  
  # Get genotype and filter
  df = lapply(vcf_l, 
               function(x) {
                 gt = geno(x)$GT %>% 
                   as.data.frame %>% 
                   rownames_to_column(var = "variant")
                 gt$filter <- VariantAnnotation::fixed(x)$FILTER
                 colnames(gt) <- c("variant","genotype",'filter')
                 return(gt)
               }) %>% 
    bind_rows(.id = "sample")

  # Add patient data
  df <- left_join(df, mito_reference_df, by = 'sample')
  df <- filter(df, complete.cases(df))
  df$varsamp <- paste(df$variant, df$patient, sep = '_')
  
  # Add frequency of variant per patient
  oc = as.data.frame(table(unlist(df$varsamp)))
  colnames(oc) <- c('varsamp','patient_freq')
  df <- left_join(df, oc, by = 'varsamp')
  df$varsamp <- NULL
  
  # Add frequency of variant overall
  oc = as.data.frame(table(unlist(df$variant)))
  colnames(oc) <- c('variant','overall_freq')
  df <- left_join(df, oc, by = 'variant')

  return(df)
}

# Load variant annotations in a tibble

get_exon_table = function(vcf) {    #Get gene annotation
  ann = info(vcf)$ANN %>%
    as.list()
  ann[elementNROWS(ann) == 0] = ""
  ann = purrr::map_chr(ann, str_c, collapse = ";")    #Filter for coding muts
  coding_f = str_detect(ann, "synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant")
  if (!sum(coding_f)){
    return(NULL)
  }
  ann_exon = ann[coding_f]    #Get info for gene table
  ann_l = str_split(ann_exon, pattern = "\\|")
  effect_type = purrr::map_chr(ann_l, 2) #Missense, nonnsense ect.
  effect = purrr::map_chr(ann_l, 3)
  gene = purrr::map_chr(ann_l, 4)
  vcf = vcf[coding_f]
  gt = geno(vcf)$GT
  tb = gt %>%
    as.data.frame() %>%
    rownames_to_column( var = 'variant') %>%
    dplyr::mutate(sample = colnames(gt[1]),
                  gene = gene,
                  effect = effect,
                  effect_type = effect_type)
  tb$sample = colnames(tb)[2]
  colnames(tb)[2] = 'genotype'
    return(tb)
}

grab_variant_annotations <- function(vcf_filter) {
  
  f <- list.files(pattern = vcf_filter, 
                  recursive = T)
  
  ff = grep(pattern = "merged", 
            x = f, 
            val = T, 
            invert = T)
  
  vcfs = lapply(ff, function(x) {
    v = readVcf(x)
    return(v)
  })
  
  names(vcfs) = paste(sapply(ff, 
                             function(x) {
                               a = str_split_fixed(x, 
                                                   pattern = "/", 
                                                   n = 2)
                               return(a[1])
                             }))
  
  annotated_variants = lapply(vcfs, get_exon_table)
  
  dfav = do.call(rbind,annotated_variants)
  dfav = dplyr::select(dfav, -genotype)
  
  return(dfav)
}

single_type_occurence <- function(grange, ref_genome) {
  df <- data.frame()
  vcf <- grange
  types <- mut_type(vcf)
  CpG = c("ACG", "CCG", "TCG", "GCG")
  column_names = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G",
                   "C>T at CpG", "C>T other")
  CT_context = 0
  CT_at_CpG = 0
  CT_at_other = 0
  CT_muts = which(types == "C>T")
  
  if (length(CT_muts) > 0) {
    CT_context = type_context(vcf[CT_muts], ref_genome)[[2]]
    CT_at_CpG = sum(!(is.na(BiocGenerics::match(CT_context,CpG))))
    CT_at_other = length(CT_muts) - CT_at_CpG
  }
  
  # Construct a table and handle missing mutation types.
  full_table = table(factor(types, levels = column_names))
  full_table["C>T at CpG"] = CT_at_CpG
  full_table["C>T other"] = CT_at_other
  df = BiocGenerics::rbind(df, full_table)
  colnames(df) = names(full_table)
  return(df)
}
