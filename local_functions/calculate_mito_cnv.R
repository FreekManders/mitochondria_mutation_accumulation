# Function: outputs DF with mean, std error, median, MAD per sample per patient


calculate_mito_cnv <- function(wd, sample_reference_df) {
  
  old_dir = setwd(wd)
  on.exit(setwd(old_dir), add = TRUE)
  
  # Find metric files
  nfiles = list.files(pattern = "WGSMetrics|wgs_metrics", recursive = T)
  mtfiles = list.files(pattern = "^metrics.txt", recursive = T)
  
  # Filter for nuclear samples in sample list
  n_sample_names = nfiles %>% 
    basename() %>% 
    str_remove(".wgs_metrics.txt|_dedup.WGSMetrics.txt")
  n_sample_f = n_sample_names %in% sample_reference_df$sample
  nfiles = nfiles[n_sample_f]
  n_sample_names = n_sample_names[n_sample_f]
  
  # Filter for mt samples in sample list
  mt_sample_names = str_split_fixed(mtfiles, pattern = "/", n = 3)[,2]
  mt_sample_f = mt_sample_names %in% sample_reference_df$sample
  mtfiles = mtfiles[mt_sample_f]
  mt_sample_names = mt_sample_names[mt_sample_f]
  
  # Read data
  ncov = lapply(nfiles, read.table, header = T, stringsAsFactors = F, nrows = 1)
  mtcov = lapply(mtfiles, read.table, header = T, stringsAsFactors = F, nrow = 1)

  
  # Extract mean, median, SD, MAD
  #Nuclear
  n_df = ncov %>% 
    bind_rows() %>% 
    dplyr::select("n_mean" = MEAN_COVERAGE, 
                  "n_sd" = SD_COVERAGE, 
                  "n_median" = MEDIAN_COVERAGE, 
                  "n_MAD" = MAD_COVERAGE) %>% 
    dplyr::mutate("sample" = n_sample_names)
  
  #MT
  mt_df = mtcov %>% 
    bind_rows() %>% 
    dplyr::select("mt_mean" = MEAN_COVERAGE, 
                  "mt_sd" = SD_COVERAGE, 
                  "mt_median" = MEDIAN_COVERAGE, 
                  "mt_MAD" = MAD_COVERAGE) %>% 
    dplyr::mutate("sample" = mt_sample_names)
  

  
  # Merge tables and compute CNV
  df <- left_join(mt_df, n_df, by = "sample")
  df <- mutate(df, cnv_mean = 2 * mt_mean / n_mean) 
  df <- mutate(df, cnv_median =  2 * mt_median / n_median)
  
  # Add sample reference table
  df <- left_join(sample_reference_df, df, by = "sample")

  # Check there are no duplicated samples
  duplicated_samples = df$sample[duplicated(df$sample)]
    
  if (length(duplicated_samples)){
    stop(paste0("The following samples are present more than once in your data: ", 
                paste0(duplicated_samples, collapse = ", ")))
  }
  
  return(df)
}
