#! /bin/bash

#SBATCH --job-name=mito_pipeline_run
#SBATCH -t 03:00:00
#SBATCH --mem=20G
#SBATCH --gres=tmpspace:20G
#SBATCH -o /hpc/pmc_vanboxtel/projects/Freek_mito/mito_pipeline_output.txt
#SBATCH -e /hpc/pmc_vanboxtel/projects/Freek_mito/mito_pipeline_errors.txt
#SBATCH --mail-user=f.m.manders@prinsesmaximacentrum.nl
#SBATCH --mail-type=ALL

#module load autossh/1.4e
#autossh -M 13306 horus -L 3306:localhost:3306 -N -f &
#autossh -M 18009 horus -R 8007:localhost:8007 -N -f &
#pid=$!
#sleep 5

# Take previous output
MAIN_LOC=$1
WDL_LOC=$2
INPUT_LOC=$3
NAME=$4

#set -o errexit
#set -o pipefail

# Get the cromwell password
source ~/.cromwell_creds.sh

# Check whether sample has been successfully run already
if ! [[ -d "$MAIN_LOC$*/final_output/$NAME" ]]
    then

        curl -X POST "https://fmanders.hpccw.op.umcutrecht.nl//api/workflows/v1" \
            --write-out \\n%{http_code} \
            --user "fmanders:${CROMWELL_PASSWORD}" \
            --header "accept: application/json" \
            --header "Content-Type: multipart/form-data" \
            --form "workflowSource=@${MAIN_LOC}${WDL_LOC}MitochondriaPipeline_JEdit.wdl" \
            --form "workflowInputs=@${MAIN_LOC}${INPUT_LOC}mitochondrial_workflow_inputs_${NAME}.json" \
            --form "workflowDependencies=@${MAIN_LOC}${WDL_LOC}modules.zip" \
            --form "workflowOptions=@${MAIN_LOC}${INPUT_LOC}mitochondrial_workflow_options_${NAME}.json"
fi
