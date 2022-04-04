get_somatic_variant_df = function(ref_genome, vcf_wd, vcf_filter, sample_reference_df, low_vaf_filter, high_vaf_filter){
    
    # Read vcfs.
    vcf_l = read_all_vcfs(wd = vcf_wd,
                          vcf_filter = vcf_filter,
                          sample_reference_df = sample_reference_df)
    
    # Filter for PASSED variants of vcfs
    passed_vcflist <- filter_passed_variants(vcf_l = vcf_l)
    
    # Create a DF with ALL called variants. This includes not passed variants,
    # which are used for somatic calling later on.
    all_variants <- grab_all_variants(vcf_l,
                                      mito_reference_df = sample_reference_df)
    
    vaf_df <- get_vaf_df(vcf_list = passed_vcflist)
    
    variant_df <- create_variant_df(vcf_list = passed_vcflist, 
                                    mito_reference_df = sample_reference_df)
    
    filtered_variants <- filter_variants(variant_df = variant_df,
                                         all_variants = all_variants,
                                         sample_reference_df = sample_reference_df)
    
    
    somatic_variant_df <- combine_mito_data(filtered_variants = filtered_variants,
                                            vaf_df = vaf_df,
                                            low_vaf_filter = low_vaf_filter,
                                            high_vaf_filter = high_vaf_filter)
    return(somatic_variant_df)
}
