#! /bin/bash

#SBATCH --job-name=mitopipeline_start
#SBATCH -t 01:00:00
#SBATCH --mem=8G
#SBATCH -o /hpc/pmc_vanboxtel/projects/Freek_mito/mitopipeline_startup_output.txt
#SBATCH -e /hpc/pmc_vanboxtel/projects/Freek_mito/mitopipeline_startup_error.txt
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL

set -o errexit

# Check existence of CONFIG file in root dir
source mitopipeline.config

##Load required modules
module load Java/1.8.0_60
module load python/3.6.1

cd $HG38_MAP_LOC

# Make sure we can look for files recursively
shopt -s globstar

FLOW=(q)
for BAMLOC in **/*.bam; do
# Remove BAM files created by structural variant calling
    if ! [[ "${BAMLOC}" =~ "gridss" || "${BAMLOC}" =~ "total_BAMS" || "${BAMLOC}" =~ "10925" || "${BAMLOC}" =~ "9379" || "${BAMLOC}" =~ "14788" || "${BAMLOC}" =~ "14436" || "${BAMLOC}" =~ "12030" || "${BAMLOC}" =~ "11498" || "${BAMLOC}" =~ "11266" || "${BAMLOC}" =~ "20904AML4D13" || "${BAMLOC}" =~ "7125DX2AML" || "${BAMLOC}" =~ "_telbam" ]]; then

# Only pick BAM files located in /mapping/ or in /BAMS/
        if [[ "${BAMLOC}" =~ "BAMS" || "${BAMLOC}" =~ "bams" || "${BAMLOC}" =~ "mapping" || "${ALL_BAMS}" == "TRUE" ]]; then
#        if [[ "${BAMLOC}" =~ "BAMS" || "${ALL_BAMS}" == "TRUE" ]]; then
            LOC=${BAMLOC%.*}/
            ARRAY=(` echo $LOC | sed "s:/: :g"`)
            BAM=${ARRAY[@]: -1:1}.bam

# Grab Sample ID from BAM file
            if ! [[ "$BAM" =~ "sorted" ]]; then
                if [[ "$BAM" =~ "_dedup" ]]; then
                    SAMPLE=`echo ${ARRAY[@]: -1:1} | sed -e "s/_dedup//"`
                    if [[ "$SAMPLE" =~ ".realigned" ]]; then
                        SAMPLE=`echo $SAMPLE | sed -e "s/.realigned//"`
                        if ! [[ " ${FLOW[*]} " =~ " ${SAMPLE} " ]]; then
                            FLOW+=($SAMPLE)

                        fi

                    else
                        SAMPLE=$SAMPLE
                        if ! [[ " ${FLOW[*]} " =~ " ${SAMPLE} " ]]; then
                            FLOW+=($SAMPLE)
                        fi
                    fi


                elif [[ "$BAM" =~ ".realigned" ]]; then
                    SAMPLE=`echo ${ARRAY[@]: -1:1} | sed -e "s/.realigned//"`
                    if ! [[ " ${FLOW[*]} " =~ " ${SAMPLE} " ]]; then
                        FLOW+=($SAMPLE)
                    fi

                else
                    SAMPLE=`echo ${ARRAY[@]: -1:1}`
                    if ! [[ " ${FLOW[*]} " =~ " ${SAMPLE} " ]]; then
                        FLOW+=($SAMPLE)
                    fi
                fi
                if [[ "${WGS_METRICS_LOC}" == "NA" ]]; then
                    if [[ "$BAM" =~ "_dedup" ]]; then # Ensures the correct wgs metric file is taken for samples that were analyzed twice.
                        QC_L=$(find -L . -maxdepth 4 -type f | egrep '*.wgs_metrics.txt' | grep -v "/._")
                        echo $SAMPLE
                        if [[ -z "$QC_L" ]]; then
                            QC_L=$(find -L . -maxdepth 4 -type f | egrep '*.wgs_metrics.txt' | grep -v "/._")
                        fi
                    else
                        QC_L=$(find -L . -maxdepth 4 -type f | egrep '*.wgs_metrics.txt' | grep -v "/._")
                        if [[ -z "$QC_L" ]]; then
                            QC_L=$(find -L . -maxdepth 4 -type f | egrep '*_dedup_WGSMetrics.txt' | grep -v "/._")
                        fi
                    fi
                else
                    QC_L=$(find -L ${WGS_METRICS_LOC} -maxdepth 1 -type f | egrep '*.wgs_metrics.txt|*_dedup_WGSMetrics.txt' | grep -v "/._")
                fi
                QC_FOUND="0"
                for QC in ${QC_L}; do
                    if [[ "$QC" =~ "$SAMPLE" && ! -d "$WORK_LOC/final_output/$SAMPLE" ]]; then
                        QC_FOUND="1"
                        cp $QC ${MAIN_LOC}/Nuclear_wgs_metrics/
                        #QCLOC=$HG38_MAP_LOC$QC
                        cd $MAIN_LOC$INPUT_LOC
    # Finish input.json with specified files
                        echo python3 finish_json_slurm.py -b $BAMLOC -n $SAMPLE -q $QC -l $HG38_MAP_LOC -wl $WORK_LOC
                        python3 finish_json_slurm.py -b $BAMLOC -n $SAMPLE -q $QC -l $HG38_MAP_LOC -wl $WORK_LOC
                        cd $HG38_MAP_LOC
                    fi
                done
                if [[ ${QC_FOUND} == "0" && ! -d "$WORK_LOC/final_output/$SAMPLE" ]]; then
                    echo "No QC file found for: $SAMPLE"
                fi
            fi
        fi
    fi
done

cd $MAIN_LOC$WDL_LOC

zip modules.zip MitochondriaPipeline_JEdit.wdl AlignmentPipeline_JEdit.wdl AlignAndCall_JEdit.wdl

cd $MAIN_LOC

echo ${FLOW[@]1}

## RUN cromwell for each specified sample ; skip first variable ( q )
for NAME in ${FLOW[@]:1}; do

    if ! [[ -d "$WORK_LOC/final_output/$NAME" ]]; then
        sbatch mito_pipeline_run2.sh $MAIN_LOC $WDL_LOC $INPUT_LOC $NAME
        sleep 5
	echo "This doesn't exist yet:"
	echo "$WORK_LOC/final_output/$NAME"
    fi
done
