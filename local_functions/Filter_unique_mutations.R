# Filter for unique mutations in list of VCFs

create_variant_df <- function(vcf_list, mito_reference_df) {
  gt_l = lapply(vcf_list, 
               function(vcf) {
                 gr <- granges(vcf)
                 
                 # Remove variants with multiple alternative alleles.
                 multi_alts_f = elementNROWS(gr$ALT) != 1
                 vcf = vcf[!multi_alts_f]
                 gr = gr[!multi_alts_f]
                 
                 # Get position, ref and alt
                 pos = start(gr)
                 ref = as.vector(gr$REF)
                 alt = as.vector(unlist(gr$ALT))
                 
                 # Get genotype
                 gt = geno(vcf)$GT[,1, drop = T]
                 
                 # Combine results
                 res = tibble("mut_pos" = pos, "ref" = ref, "alt" = alt, "genotype" = gt, "variant" = names(vcf))
                 return(res)
               })
  
  # Combine samples
  geno_df = bind_rows(gt_l, .id = "sample")
  
  # Add patient data
  geno_df <- left_join(geno_df, mito_reference_df, by = 'sample')
  
  geno_df = split(geno_df, geno_df$patient)
  geno_df = do.call(rbind, geno_df)
  
  return(geno_df)
}

count_vars <- function(patient_df, mito_reference_df) {
  
  # Count number of occurrences of variant in patient
  spm = lapply(patient_df, function(x) {
    xx = unique(x)
    oc = as.data.frame(table(unlist(xx$variant)))
    colnames(oc) <- c('variant','freq')
    
    # Connect number of occurrences to variant
    a = left_join(xx, oc, by = "variant")
    return(a)
  })
  
  # Add them back together
  dfum = do.call(rbind,spm)
  oc = as.data.frame(table(unlist(mito_reference_df$patient)))
  colnames(oc) <- c('patient','max_freq')
  dfum = left_join(dfum, oc, by = 'patient')
  dfum = filter(dfum, dfum$genotype %in% c('0/1','0|1'))
  return(dfum)
}


filter_variants <- function(variant_df, all_variants, sample_reference_df) {
  
  # Combine variant_df and all_variants
  vdf <- left_join(variant_df, all_variants)
  
  # Determine nr. samples per patient/donor This is the maximum frequency of that variant
  oc <- as.data.frame(table(unlist(sample_reference_df$patient)))
  colnames(oc) <- c('patient','max_freq')
  vdf <- left_join(vdf, oc)
  
  # Filter variants
  # Filter for correct genotype
  vdf <- filter(vdf, genotype %in% c('0/1','0|1')) %>%
    # Filter for occurrence in patient
    filter(max_freq == 1 | patient_freq < max_freq) %>%
    # Filter for global occurrence. A variant can't occur in more than one patient
    filter(overall_freq <= patient_freq)
  
  # Remove variants that are shared and for which there is evidence in a bulk sample
  bulk_shared_variants = dplyr::filter(all_variants, bulk == "bulk" & patient_freq > 1)
  vdf = anti_join(vdf, bulk_shared_variants, by = c("patient", "genotype", "variant"))
  
  return(vdf)
}
