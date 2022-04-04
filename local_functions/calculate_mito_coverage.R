# Calculate Mitochondrial coverage


# supply list of of patients and the dir that contains the tabular coverage data
calculate_mt_coverage <- function(wd, mito_reference_df) {
  old_dir = setwd(wd)
  on.exit(setwd(old_dir), add = TRUE)
  
  # Find base coverage files
  mtfiles = list.files(pattern = "per_base_coverage.tsv", all.files = T, recursive = T)
  
  # Remove files not matching the sample reference table
  mt_sample_names = str_split_fixed(mtfiles, pattern = "/", n = 3)[,2]
  mt_sample_f = mt_sample_names %in% mito_reference_df$sample
  mtfiles = mtfiles[mt_sample_f]
  mt_sample_names = mt_sample_names[mt_sample_f]

  
  # Read data
  df_l = lapply(mtfiles, read.table, header = T)
  df = bind_rows(df_l, .id = "sample")
  
  
  # Add sample names
  df$sample = rep(mt_sample_names, elementNROWS(df_l))
  
  # Remove column and change names
  df$target = NULL
  colnames(df) = c("sample",'chrom','position','coverage')
  
  # Add sample reference file
  rdat <- left_join(df, mito_reference_df, by = "sample")
  
  return(rdat)
}
