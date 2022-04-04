combine_mito_data <- function(filtered_variants, 
                              vaf_df, 
                              low_vaf_filter, 
                              high_vaf_filter) {
  
  # Join informational DFs
  cdf = left_join(filtered_variants, vaf_df)
  df = filter(cdf, 
                    cdf$VAF >= low_vaf_filter & cdf$VAF <= high_vaf_filter)
  
  return(df)
}
