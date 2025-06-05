#!/usr/bin/env bash

SCRIPT_PATH=$(dirname "$(realpath "$0")") #Get script's path to source globalFunctions
source ${SCRIPT_PATH}/globalFunctions.sh #Import custom functions eg. createOutDir

while getopts "i:o:r:" opt; do
    case $opt in
        i) merged_vcf="$OPTARG" ;;  #/path/to/unfiltered/vcf/file
        o) output_dir="$OPTARG" ;;  #/path/to/output/folder - can be new or existing
        r) reference="$OPTARG" ;;   #/path/to/reference.fasta
        *)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i <path/to/VCF> -o <path/to/output> -r <reference.fasta>"
            exit 1
            ;;
    esac
done

checkInputArgs merged_vcf output_dir reference

if [ -f ${merged_vcf} ]; then
    sample_name=$(basename ${merged_vcf} .vcf.gz )
    echo "Found ${sample_name}. Proceeding with hard filtering"
fi

log_file="${sample_name}_variantHardFilter.log"

createOutDir ${output_dir} #Create output directory

processLog ${output_dir}/${log_file} && start_time=$(logrunTime)

# bcftools is needed to reorder variants by position 
checkProgram gatk bcftools parallel

# Define pipeline functions
sort_VCF(){
    local input_vcf=$1 #${sample_name}.vcf.gz
    local sortedVCF=$2 #${sample_name}_${variant}.vcf.gz
    bcftools sort ${input_vcf} -o ${sortedVCF}
    exit_value=$?
    checkExitStatus "bcftools sort" $?
}

reIndex_VCF(){
    local vcf=$1
    echo "[INFO] Re-indexing ${vcf}"
    gatk IndexFeatureFile   \
        -I ${vcf}
    exit_value=$?
    checkExitStatus "IndexVCF" $?
}

run_SelectVariants(){
    local variant_type=$1 #SNP or INDEL
    echo "[INFO] Separating ${variant_type}s from main VCF"
    gatk SelectVariants \
        -R ${reference} \
        -V ${sorted_vcf} \
        -select-type ${variant_type} \
        -O ${output_dir}/${sample_name}_${variant_type}.vcf.gz
    checkExitStatus "SelectVariants ${variant_type}" $?
}
#############################################
# Adjust parameters for filtering below     #
# Default GATK filtering parameters are set #
#############################################

run_VariantFiltration(){
    local variant_type=$1
    #Output from SelectVariants
    local input=${output_dir}/${sample_name}_${variant_type}.vcf.gz
    echo "[INFO] Running VariantFiltration for ${variant_type}"
    if [[ ${variant_type} == "SNP" ]]; then
        gatk VariantFiltration  \
            -V ${input} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O ${output_dir}/${sample_name}_SNP_filtered.vcf.gz
        checkExitStatus "VariantFiltration SNP" $? 
    elif [[ ${variant_type} == "INDEL" ]]; then
        gatk VariantFiltration \
            -V ${input} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O ${output_dir}/${sample_name}_INDEL_filtered.vcf.gz
        checkExitStatus "VariantFiltration INDEL" $?
    fi
}

# This is where the pipeline runs based on the created functions above #
sorted_vcf=${output_dir}/${sample_name}_sorted.vcf.gz
sample_StartRunTime=$(logrunTime)
sort_VCF ${merged_vcf} ${sorted_vcf}
reIndex_VCF ${sorted_vcf}

# Export functions and variables for parallel 
export -f run_SelectVariants
export -f run_VariantFiltration
export -f reIndex_VCF
export -f checkExitStatus
export sample_name output_dir sorted_vcf reference

# Run the filtering process for SNP and INDEL in parallel #
parallel --jobs 2 --halt soon,fail=1 '
    variant_type={1}
    run_SelectVariants $variant_type
    run_VariantFiltration $variant_type
    reIndex_VCF ${output_dir}/${sample_name}_${variant_type}_filtered.vcf.gz
' ::: SNP INDEL

# Merge SNPs and INDELs to one vcf file
gatk MergeVcfs \
    I=${output_dir}/${sample_name}_SNP_filtered.vcf.gz \
    I=${output_dir}/${sample_name}_INDEL_filtered.vcf.gz \
    O=${output_dir}/${sample_name}_filtered.vcf.gz
checkExitStatus "Merge VCFs" $?
sample_endRunTime=$(logrunTime)

echo "Removing intermediate files"
rm -f ${output_dir}/${sample_name}_{SNP,INDEL}.vcf.{gz,gz.tbi} \
    ${output_dir}/${sample_name}_sorted.vcf.{gz,gz.tbi}

endLog ${sample_StartRunTime} ${sample_endRunTime} ${sample_name}
checkScriptExit $?
