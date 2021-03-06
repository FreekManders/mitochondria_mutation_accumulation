import argparse
import os

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-b','--bam',type=str, help='Input BAM file', required = True)
parser.add_argument('-n','--name',type=str, help='Input sample name', required = True)
parser.add_argument('-q','--wgsmetrics',type=str, help='Input location of WGSMetrics file', required = True)
parser.add_argument('-l','--location',type=str, help='Input location of \'root\' folder of BAM files', required = True)
parser.add_argument('-wl','--worklocation',type=str, help='Input work folder for analysis', required = True)
parser.add_argument('-m','--metrics_location',type=str, help='Location of wgs_metrics files', required = True)
args = parser.parse_args()

## Declare location of file which contains autosomal coverage
if (args.metrics_location == "NA"):
    cov_loc = args.location + args.wgsmetrics
else:
    cov_loc = args.wgsmetrics

index_loc = "/hpc/pmc_vanboxtel/projects/Markus_AML_ALL_tracing/hg38/Resubmission/IAP/PMC21586/bams/"
#index_loc = args.location

## Create name of index file
index_file = index_loc + args.bam + '.bai'
if not os.path.isfile(index_file):
    index_file = index_loc + args.bam[:-4] + '.bai'
if not os.path.isfile(index_file):
    index_file = index_loc + args.bam[:-4] + '_dedup.bai' #Adds a second '_dedup'
if not os.path.isfile(index_file):
    index_file = index_loc + os.path.basename(args.bam)[:-4] + ".bai"


autosomal_cov = 30 ## Default value

with open('mitochondrial_workflow_inputs_unfinished.json', mode='r') as datafile:
    input_lines = datafile.readlines()
    datafile.close()

    with open(cov_loc, mode='r') as covfile:

        cov_list = []
        for line in covfile:
            line = line.rstrip()
            splitline = line.split('\t')
            cov_list.append(splitline)
        ## With Van Boxtel IAP, this is the correct value for autosomal coverage
        autosomal_cov = cov_list[7][1]
        covfile.close()

        ## Lines to write in new JSON
        bam_loc = '   \"MitochondriaPipeline.wgs_aligned_input_bam_or_cram\" : ' + '\"' + args.location + args.bam + '\",'
        bai_loc = '   \"MitochondriaPipeline.wgs_aligned_input_bam_or_cram_index\" : ' + '\"' + index_file + '\",'
        name = '   \"MitochondriaPipeline.sample_name\" : ' + '\"' + args.name + '\",'
        cov = '   \"MitochondriaPipeline.autosomal_coverage\" : ' + '\"' + str(autosomal_cov) + '\"'
        input = [bam_loc, bai_loc, name,cov]

        ## Create new json files
        inputlogname = 'mitochondrial_workflow_inputs_' + args.name + '.json'
        optionlogname = 'mitochondrial_workflow_options_' + args.name + '.json'

        ## Create new input file for sample specific json
        inputlog = open(inputlogname, mode="w")
        for line in input_lines:
            inputlog.write(line)
        for i in input:
            inputlog.write(i)
            inputlog.write('\n')
        inputlog.write('}')
        inputlog.close()

        ## Create new options file for sample specific json
        optionlog = open(optionlogname, mode="w")
        optionlog.write('{')
        optionlog.write('\n')
        optionlog.write('   \"final_workflow_outputs_dir\": ' + '\"' + args.worklocation + 'final_output/' + args.name + '/\"')
        optionlog.write('\n')
        optionlog.write('}')
        optionlog.close()
