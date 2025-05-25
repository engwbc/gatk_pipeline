#!/usr/bin/env bash

SCRIPT_PATH=$(dirname "$(realpath "$0")") #Get script's path to source globalFunctions
source ${SCRIPT_PATH}/globalFunctions.sh #Import custom functions eg. createOutDir
exit_value= # Assign a variable for exit code $?

while getopts "i:o:r:H:O:T:M:s:" opt; do
    case $opt in
        i) merged_VCF="$OPTARG" ;;  #/path/to/vcf/file
        o) output_dir="$OPTARG" ;;  #/path/to/output/folder - can be new or existing
        r) reference="$OPTARG" ;;   #/path/to/reference.fasta
        H) HAPMAP="$OPTARG" ;;      #/path/to/hapmap.hg38.vcf.gz
        O) OMNI="$OPTARG" ;;        #/path/to/1000G_omni2.5.hg38.vcf.gz
        T) thousandGenomes="$OPTARG" ;; #/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        M) MILLS="$OPTARG" ;;       #/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        s) dbSNP="$OPTARG" ;;        #/path/to/dbsnp138.vcf.gz
        *) 
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i <merged_VCF> -o <output_dir> -r <reference.fasta> -H <hapmap.hg38> -O <1000G_omni2.5.hg38.vcf.gz> -T <1000G_phase1.vcf.gz> -M <Mills_and_1000G.vcf.gz> -s <dbSNP.vcf.gz>" >&2
            exit 1
            ;;
    esac
done

checkInputArgs merged_VCF output_dir reference

log_file="variantQC.log"

createOutDir ${output_dir} #Create output directory

processLog ${output_dir}/${log_file} && start_time=$(logrunTime)

# bcftools is needed to reorder variants by position 
checkProgram gatk bcftools

printf "Input directory: ${merged_VCF}\nOutput directory: ${output_dir}\n"

if [ -f ${merged_VCF} ]; then
    sample_name=$(basename ${merged_VCF} .vcf.gz )
    echo "Found ${sample_name}"
    echo "Proceeding with VariantRecalibrator"
    printf "\n[INFO] For VariantRecalibrator to work well, the data should at least\
 be a single whole-genome or 30 exomes (humans). Anything smaller may cause inaccurate indel recalibration.\n" 
else
    echo "[ERROR] Could not find VCF file"
    end_time=$(logrunTime) && endLog "$start_time" "$end_time"
    exit 2
fi

## Reorder reads and create a tabix index file for the merged VCF
sorted_vcf=${output_dir}/${sample_name}_sorted.vcf.gz
bcftools sort ${merged_VCF} -o ${sorted_vcf}
exit_value=$?
checkExitStatus "bcftools sort" $exit_value

gatk IndexFeatureFile   \
    -I ${sorted_vcf}
exit_value=$?
checkExitStatus "IndexVCF" $exit_value

##                                  ##
## Run VariantRecalibrator for SNPs ##
##                                  ##

# Make sure to change the path for each of the training/truth resource sets
#Uses hapmap as the highest confidence resource and thus representative of true sites (truth=true)
#1000G Omni and Phase1 are used as training sets for the recalibration model (training=true,truth=false)
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
    -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -R ${reference} \
    -V ${sorted_vcf}  \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0   \
    ${HAPMAP}    \
    --resource:omni,known=false,training=true,truth=false,prior=12.0    \
    ${OMNI} \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
    ${thousandGenomes} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  \
    -mode SNP -O ${output_dir}/${sample_name}_SNP.recal \
    --tranches-file ${output_dir}/${sample_name}_SNP.tranches \
    --rscript-file ${output_dir}/${sample_name}_SNP.plots.R
exit_value=$?
checkExitStatus "SNP Recalibration" $exit_value

##                                              ##
##      Run VariantRecalibrator for Indels      ##
##                                              ##

# Make sure to change the path for each of the training/truth resource sets
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
    -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -R ${reference} \
    -V ${sorted_vcf}  \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 \
    ${MILLS} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
    ${dbSNP} \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL -O ${output_dir}/${sample_name}_indels.recal \
    --tranches-file ${output_dir}/${sample_name}_indels.tranches \
    --rscript-file ${output_dir}/${sample_name}_indels.plots.R
exit_value=$?
checkExitStatus "Indel Recalibration" $exit_value

##                              ##
##      ApplyVQSR for SNPs      ##
##                              ##

# Define separate variables for SNP and SNP+Indel recalibrated VCFs
snp_recalibratedVCF=${output_dir}/${sample_name}.SNP.recalibrated_99.9.vcf.gz
snp_indel_recalibratedVCF=${output_dir}/${sample_name}.filtered.vcf.gz

#From VariantRecalibrator, we select for variants in the top three tranches (100, 99.95 and 99.9)
#This corresponds to allowing a 0.1% rate of false positives
# This is preliminary and adjustable based on what works best for your sample
gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V ${sorted_vcf}    \
    --recal-file ${output_dir}/${sample_name}_SNP.recal \
    --mode SNP  \
    --tranches-file ${output_dir}/${sample_name}_SNP.tranches \
    --truth-sensitivity-filter-level 99.9   \
    --create-output-variant-index true  \
    -O ${snp_recalibratedVCF}
exit_value=$?
checkExitStatus "ApplyVQSR SNP" $exit_value

## ApplyVQSR doesn't like it if the VCF is not indexed, so this is repeated for safety measures
gatk IndexFeatureFile   \
    -I ${snp_recalibratedVCF}
exit_value=$?
checkExitStatus "IndexRecalVCF SNP" $exit_value

##                                  ##
##      ApplyVQSR for Indels        ##
##                                  ##
gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V ${snp_recalibratedVCF}   \
    --recal-file ${output_dir}/${sample_name}_indels.recal   \
    --mode INDEL    \
    --tranches-file ${output_dir}/${sample_name}_indels.tranches \
    --truth-sensitivity-filter-level 99.9   \
    --create-output-variant-index true  \
    -O ${snp_indel_recalibratedVCF}
checkExitStatus "ApplyVQSR INDEL" $exit_value

gatk IndexFeatureFile   \
    -I ${snp_indel_recalibratedVCF}
exit_value=$?
checkExitStatus "IndexRecalVCF INDEL" $exit_value

checkScriptExit $exit_value
