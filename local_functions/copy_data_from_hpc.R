copy_data_from_hpc = function(source_dir, dest_dir, sample_folder, sample_reference_df){
    
    # Go to source dir
    old_dir = setwd(source_dir)
    on.exit(setwd(old_dir), add = TRUE)
    
    # Copy mito coverage data
    samples_final_output = list.files(file.path(sample_folder, "final_output"), full.names = TRUE)
    
    # Remove samples not in sample list
    samples_final_output = samples_final_output[basename(samples_final_output) %in% sample_reference_df$sample]
    purrr::map(samples_final_output, copy_mito_metrics, dest_dir)
    
    # Copy WGS coverage data
    copy_wgs_metrics(dest_dir, samples_final_output)
    
    # Copy annotated vcfs
    copy_vcfs(dest_dir, samples_final_output)
}

copy_mito_metrics = function(sample_final_output, dest_dir){
    
    # Determine destination directory
    sample = basename(sample_final_output)
    dest_dir = file.path(dest_dir, "output", "mito_coverage", sample)
    if (!dir.exists(dest_dir)){
        dir.create(dest_dir)
    }
    
    # Create list of files to copy
    files = list.files(sample_final_output, full.names = TRUE)
    files = str_subset(files, ".\\.bam$|.\\.vcf$", negate = TRUE)
    files_to = file.path(dest_dir, basename(files))
    
    # Copy files
    file.copy(files, files_to, overwrite = FALSE)
    return(0)
}

copy_wgs_metrics = function(dest_dir, samples_final_output){

    # Find files
    wgs_metrics = list.files("Nuclear_wgs_metrics", full.names = TRUE)
    
    # Remove samples not in sample list
    sample_names = wgs_metrics %>% 
        basename() %>% 
        str_remove("_dedup_WGSMetrics.txt|.wgs_metrics.txt")
    wgs_metrics = wgs_metrics[sample_names %in% basename(samples_final_output)]
    
    # Determine destination
    wgs_metrics_dest = file.path(dest_dir, "output", "nuclear_coverage", basename(wgs_metrics))
    
    # Copy files
    file.copy(wgs_metrics, wgs_metrics_dest, overwrite = FALSE)
    return(0)
}

copy_vcfs = function(dest_dir, samples_final_output){
    
    # Find files
    vcf_fnames = list.files(file.path("annotated_vcfs", "VCFS", "VCF"), full.names = TRUE)
    
    # Remove samples not in sample list
    samples = vcf_fnames %>% 
        basename() %>% 
        str_remove("\\.filtered.*")
    sample_f = samples %in% basename(samples_final_output)
    vcf_fnames = vcf_fnames[sample_f]
    samples = samples[sample_f]
    
    # Determine destination directories and create them
    dest_dirs = file.path(dest_dir, "annotated_vcfs", samples)
    dest_dirs_create = unique(dest_dirs[!dir.exists(dest_dirs)])
    if (length(dest_dirs_create)){
        purrr::map(dest_dirs_create, dir.create)
    }
    
    # Determine destination files
    vcf_fnames_dest = file.path(dest_dirs, basename(vcf_fnames))
    
    # Copy files
    file.copy(vcf_fnames, vcf_fnames_dest, overwrite = FALSE)
    return(0)
}
