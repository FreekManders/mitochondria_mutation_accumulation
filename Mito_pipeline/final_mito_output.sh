#!/bin/bash

#SBATCH --job-name=finalise_output_mito_pipeline

source mitopipeline.config


cd $WORK_LOC
for FOLDER in final_output/*; do
    echo $WORK_LOC$FOLDER
    cd $WORK_LOC$FOLDER

    SAMPLE_NAME=$(basename ${FOLDER})
    GREP_STRING="theoretical_sensitivity.txt$|metrics.txt$|"${SAMPLE_NAME}"_dedup.realigned.contamination.txt$|"${SAMPLE_NAME}".vcf.idx$|"${SAMPLE_NAME}".vcf$|."${SAMPLE_NAME}"_dedup.realigned.bam$|"${SAMPLE_NAME}"_dedup.metrics$|"${SAMPLE_NAME}"_dedup.realigned.bai$|per_base_coverage.tsv$|"${SAMPLE_NAME}"_dedup.bai$|"${SAMPLE_NAME}"_dedup.bam$"
    find . -mindepth 2 -type f -print0 | egrep -zZ ${GREP_STRING} | xargs -0 -I{} mv {} .
    #find . -mindepth 2 -type f -print0 | egrep -zZ "theoretical_sensitivity.txt$|metrics.txt$|.*_dedup.realigned.contamination.txt$|.*.vcf.idx$|.*.vcf$|.*_dedup.realigned.bam$|.*_dedup.metrics$|.*_dedup.realigned.bai$|per_base_coverage.tsv$|.*_dedup.bai$|.*_dedup.bam$" | xargs -0 -I{} mv {} .

    cd $WORK_LOC
done
