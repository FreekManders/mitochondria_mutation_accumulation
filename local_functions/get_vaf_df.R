# Function: Get VAF in df


get_vaf_df <- function(vcf_list) {
  
  # Adds point mutation context
  dfvs = lapply(vcf_list, function(x) { 
    ad = geno(x)$AD[,1,drop=T]

    ref_ad = purrr::map_dbl(ad, 1)
    alt_ad = purrr::map_dbl(ad, 2)
    vaf = alt_ad / (ref_ad + alt_ad)
    vaf[is.na(vaf)] = 0 #For sites with 0 reads

    vaf_df = data.frame("VAF" = vaf, "variant" = names(x), "ref_ad" = ref_ad, "alt_ad" = alt_ad)
    return(vaf_df)
    })
  
  
  dfv = bind_rows(dfvs, .id = "sample") %>% 
    dplyr::select(VAF, sample, variant, ref_ad, alt_ad)
  rownames(dfv) = seq_len(nrow(dfv))

  return(dfv)
}
