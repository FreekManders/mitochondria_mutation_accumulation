#!/bin/bash

#SBATCH --job-name=annotate_chrM_to_MT
#SBATCH -o vcf_convert_log.out
#SBATCH -e vcf_convert_errlog.out

source mitopipeline.config

cd $MAIN_LOC

# Make sure we can look for files recursively
shopt -s globstar

for VCF in */*/*/*.vcf; do
    LOC=${VCF%.*}/
    ARRAY=(` echo $LOC | sed "s:/: :g"`)
    NEWVCF=${ARRAY[@]: -1:1}.vcf
    SAMPLE=${ARRAY[@]: -1:1}

# Convert chrM to MT in all VCFs
    if ! [[ -d "$MAIN_LOC/annotated_vcfs/$SAMPLE" ]]
    then
        awk '{gsub(/^chrM/,"MT"); print}' $MAIN_LOC$VCF > ${MAIN_LOC}annotated_input/temp_$NEWVCF
        awk '{gsub(/##contig=<ID=chrM,length=16569/,"##contig=<ID=MT,length=16569,assembly=Homo_sapiens_assembly38.fasta"); print}' ${MAIN_LOC}annotated_input/temp_$NEWVCF > ${MAIN_LOC}annotated_input/$NEWVCF
        rm ${MAIN_LOC}annotated_input/temp_$NEWVCF
# Add index
        /hpc/local/CentOS7/pmc_vanboxtel/bin/IGVTools/igvtools index ${MAIN_LOC}annotated_input/$NEWVCF


    fi
done

# Run annotation for all files
python3 /hpc/pmc_vanboxtel/projects/Freek_mito/create_annotate_vcf_config.py -v ${MAIN_LOC}/annotated_input/ -o $MAIN_LOC/annotated_vcfs/
bash start_NF-IAP.sh
