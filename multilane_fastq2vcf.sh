#!/usr/bin/env bash
SCRIPT_PATH=$(dirname "$(realpath "$0")") #Get script's path to source globalFunctions
source ${SCRIPT_PATH}/globalFunctions.sh #Import custom functions eg. createOutDir
exit_value= #Initialise and assign a variable for exit code $? 

while getopts "i:o:r:s:L:h" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;    #/path/to/fastq/folder
        o) output_dir="$OPTARG" ;;   #/path/to/output/folder
        r) reference="$OPTARG" ;;    #/path/to/reference/ref.fasta
        s) dbsnp="$OPTARG" ;;        #/path/to/knownSNPs.vcf.gz - for BaseRecalibration
        L) intervals="$OPTARG";;     #/path/to/variant_call_interval.file; in the future if scattered intervals exist, create a new option to skip SplitIntervals
        h) 
            printUsage "$@"
            exit 0 ;;
        *) 
            echo "Invalid option: -$OPTARG" >&2 #Handle invalid options
            printUsage "$@"
            exit 1 #Fail
            ;;
    esac
done

#Check that all arguments are given
checkInputArgs

log_file="parallelRun.log" #for record-keeping and debugging

#Check and create output directory if missing
createOutDir $output_dir

#Inititate logging
processLog ${output_dir}/${log_file} && start_time=$(logrunTime)
num_threads=$(nproc)
echo "[INFO] Number of CPU threads available: $num_threads"

#Check for installed program dependencies using the checkProgram function
checkProgram bwa gatk samtools fastp java parallel

#For logging purposes
printf "\nInput directory: ${input_dir}\nReference: ${reference}\nOutput directory: ${output_dir}\n"

split_int_dir=${output_dir}/tmp
SCATTER_COUNT=24
createTmpDir $split_int_dir
echo "Running SplitIntervals..."
gatk SplitIntervals \
    -R ${reference} \
    -L ${intervals} \
    --scatter-count ${SCATTER_COUNT} \
    -O ${split_int_dir}
exit_value=$? #get gatk exit code
checkExitStatus "SplitIntervals" $exit_value #if exit code = 0 - success, else end script
scattered_intervals=()

# Loop through each .interval_list file in the directory
for file in ${split_int_dir}/*.interval_list; do #can also be .bed, .vcf or .list
    # Check if the file exists to handle cases where no matching files are found
    if [[ -f "$file" ]]; then
        scattered_intervals+=("$file") #Append matched interval file to array
    fi
done
#Sanity check to confirm that interval lists within the directed path is equal to SCATTER_COUNT
if [ ${#scattered_intervals[@]} -eq ${SCATTER_COUNT} ]; then
    echo "SUCCESS: Split into ${#scattered_intervals[@]} intervals"
else #Direct the user to check split_int_dir in case of user error
    echo "ERROR: Expected ${SCATTER_COUNT} interval files but found ${#scattered_intervals[@]}, please check $split_int_dir"
    exit 1
fi 

#Index reference
checkRefIndex $reference
#if reference hasn't been indexed in a while or if not frequently used, safer to reindex

## Retrieve fastQ files ##
echo "[INFO] Searching for samples in ${input_dir}."
echo "[INFO] find command is set to -maxdepth 1. To avoid errors, provide the full path to the FOLDER containing the fastQ file(s)."

#for converting samples where more than 1 lane was used
# in non-interleaved PE files (forward/reverse reads are in separate fastq files)
fastq=$(find ${input_dir} -maxdepth 1 -type f \( -name "*L001_R1*.fq.gz" -o -name "*L001_R1*.fastq.gz" \))
sample_name=$(basename $fastq | sed 's/_L.*//')
echo "Found: ${sample_name}"

fastq1=${output_dir}/${sample_name}_concat_1.fq.gz
fastq2=${output_dir}/${sample_name}_concat_2.fq.gz

# as fastq files contain unsorted raw reads (no headers or metadata)
# we can concatenate the forward and reverse reads from different lanes
lanes=$(find ${input_dir} -type f -name ${sample_name}*_L*_R1_*.fq.gz -o -name ${sample_name}*_L*_R1_*.fastq.gz | sed -E 's/.*_L([0-9]{3})_.*/\1/' | sort -u)

#prepare concatenated files
> ${fastq1} 
> ${fastq2} 
# iterate through each lane and assign reads to either fowrard (read1) or reverse (read2)
for lane in ${lanes}; do
    read1=$(find ${input_dir} -type f -name "${sample_name}*_L${lane}_R1_*.fq.gz" -o -name "${sample_name}*_L${lane}_R1_*.fastq.gz")
    read2=$(find ${input_dir} -type f -name "${sample_name}*_L${lane}_R2_*.fq.gz" -o -name "${sample_name}*_L${lane}_R2_*.fastq.gz")
    # concatenate the reads per lane to a single fastq file
    if [ -n ${read1} ] && [ -n ${read2} ]; then
        echo "Adding lane ${lane}: ${read1} and ${read2}"
        cat $read1 >> ${fastq1} #assign to read1
        cat $read2 >> ${fastq2} #assign to read2
    else
        if [ -z $read1 ]; then
            echo "ERROR: Missing forward read for lane ${lane} in ${sample_name}"
            exit 1
        fi
        if [ -z $read2 ]; then
            echo "ERROR: Missing reverse read for lane ${lane} in ${sample_name}"
            exit 1
        fi
    fi  
done 
# to check that concatenation worked
#https://www.biostars.org/p/81924/#9520618 - md5 sum checked

QC_DIR="${output_dir}/fastp_out"
createOutDir $QC_DIR

######################################
### Step 1. Sequence QC with fastp ###
######################################

echo "Running fastp for adapter trimming and read QC..."
fastp \
-i ${fastq1} \
-I ${fastq2} \
-o ${QC_DIR}/${sample_name}_qc_1.fq.gz \
-O ${QC_DIR}/${sample_name}_qc_2.fq.gz \
-h ${QC_DIR}/${sample_name}_qc.html \
-j ${QC_DIR}/${sample_name}_qc.json \
--thread 4

exit_value=$?
checkExitStatus "fastp" $exit_value
echo "[INFO] Please check fastp output with multiqc!"

fastq1_qc=${QC_DIR}/${sample_name}_qc_1.fq.gz
fastq2_qc=${QC_DIR}/${sample_name}_qc_2.fq.gz

################################################
### Step 2. Alignment and positional sorting ###
################################################

echo "Starting BWA-MEM alignment and sorting for ${sample_name}"
printf "\n[INFO] Expected bwa mem stdout: [M::mem_pestat]. \
If [E::bns_fetch_seq] is shown, your reference index might be corrupt.\
Run bwa index and try again.\n" 
bwa mem -M -t 8 -v 2 \
    ${reference} \
    ${fastq1_qc} \
    ${fastq2_qc} \
    -R "@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPL:ILLUMINA" \
| samtools sort -@ 8 -o ${output_dir}/${sample_name}_sorted.bam 
#CHANGE PLATFORM NAME IF NOT OBTAINED FROM ILLUMINA SEQUENCERS 
exit_value=$?
checkExitStatus "Alignment & Pos.Sorting" $exit_value  

echo "Starting GATK MarkDuplicates for $sample_name"
sorted_bam=${output_dir}/${sample_name}_sorted.bam

##############################
### Step 3. MarkDuplicates ###
##############################

gatk MarkDuplicates \
    -I ${sorted_bam} \
    -O ${output_dir}/${sample_name}_sorted_mkd.bam \
    -M ${output_dir}/${sample_name}_marked_dup_metrics.txt #Metrics report file
exit_value=$?
checkExitStatus "MarkDuplicates" $exit_value

#Re-index marked samples
samtools index ${output_dir}/${sample_name}_sorted_mkd.bam 
checkExitStatus "MarkDuplicates Index" $exit_value

sortedMkdup_bam=${output_dir}/${sample_name}_sorted_mkd.bam

##################################
### Step 4. Base Recalibration ###
##################################

echo "Starting GATK BaseRecalibrator for $sample_name"
gatk BaseRecalibrator \
    -I ${sortedMkdup_bam} \
    -R ${reference} \
    --known-sites ${dbsnp} \
    -O ${output_dir}/${sample_name}_recal_data.table

exit_value=$?
checkExitStatus "BaseRecalibrator" $exit_value

echo "Starting GATK ApplyBQSR for $sample_name"

bqsr_result=${output_dir}/bqsr
createOutDir $bqsr_result

#########################################################
# Step 5. Apply Base Quality Score Recalibration (BQSR) #
#########################################################

gatk ApplyBQSR \
    -I ${sortedMkdup_bam} \
    -R ${reference} \
    --bqsr-recal-file ${output_dir}/${sample_name}_recal_data.table \
    -O ${bqsr_result}/${sample_name}_sorted_mkd_bqsr.bam
exit_value=$?
checkExitStatus "ApplyBQSR" $exit_value
bqsr_bam=${bqsr_result}/${sample_name}_sorted_mkd_bqsr.bam

samtools index ${bqsr_bam}
exit_value=$?
checkExitStatus "Index BQSR" $exit_value

haplo_outputfolder=${output_dir}/haplotypecaller
createOutDir $haplo_outputfolder

# Setting up variables to export for GNU parallel
haplo_outputfile=${haplo_outputfolder}/${sample_name}
export bqsr_bam reference haplo_outputfile

################################
# Step 6. GATK HaplotypeCaller #
################################

# Define command for GNU Parallel
run_HaplotypeCaller(){
    i=$1
    interval_id=$(basename "$i" | cut -d'-' -f1) #get numeric ID from interval file
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        -I ${bqsr_bam} \
        -R ${reference} \
        -ERC GVCF \
        -L ${i} \
        -O ${haplo_outputfile}_${interval_id}.g.vcf.gz  \
        --native-pair-hmm-threads 2
}
export -f run_HaplotypeCaller
parallel --eta --jobs 3 run_HaplotypeCaller ::: "${scattered_intervals[@]}"    
#By setting 4 jobs, this should be 4 intervals per task, using 8 cores in total  
exit_value=$?
checkExitStatus "HaplotypeCaller" $exit_value

## Create sample name map file for GenomicsDBImport to merge split g.vcf files##
echo "Found `find ${haplo_outputfolder} -type f -name "*\.g.vcf.gz" | wc -l` g.vcf.gz files" #find in haplotypecaller folder
echo "Creating sample list for GenomicsDBImport"
> ${output_dir}/${sample_name}.map #for debugging purposes the .map file is cleared if it exists from a previous run for whatever reason
for file in ${haplo_outputfolder}/*.g.vcf.gz; do
    interval_id=$(basename $file | sed -E 's/.*_([0-9]{4})\.g\.vcf\.gz/\1/') #extract interval id
    echo $interval_id
    map_sample_name=${sample_name}_${interval_id}
    echo -e "${map_sample_name}\t${file}" >> ${output_dir}/${sample_name}.map
done

if [ -f "${output_dir}/${sample_name}.map" ]; then
    echo "[INFO] Finished creating ${sample_name}.map file"
else
    echo "ERROR: Could not generate .map file for GenomicsDBImport!"
    exit 1
fi

export output_dir sample_name #REMEMBER TO EXPORT VARIABLES FOR PARALLEL TO RECOGNISE!
run_GenomicsDBImport(){
    i=$1
    interval_id=$(basename "$i" | cut -d'-' -f1)
    gatk --java-options "-Xmx2g -Xms2g -XX:ParallelGCThreads=2" GenomicsDBImport \
        --genomicsdb-workspace-path ${output_dir}/tmp/${sample_name}_${interval_id} \
        --sample-name-map ${output_dir}/${sample_name}.map \
        --tmp-dir ${output_dir}/tmp \
        --merge-input-intervals \
        --bypass-feature-reader \
        -L ${i}
}
export -f run_GenomicsDBImport 
parallel --eta --jobs 3 run_GenomicsDBImport ::: ${scattered_intervals[@]}
exit_value=$?
checkExitStatus "GenomicsDBImport" $exit_value

###########################################
#   Step 7. GenotypeGVCFs (gvcf to vcf)   #
###########################################

echo "Performing joint genotyping for $sample_name from HaplotypeCaller output"
echo "[INFO] For cohort studies (e.g. families), use CombineGCVFs prior to this step to compile multiple!"
export reference output_dir sample_name

run_GenotypeGVCFs(){
    i=$1
    interval_id=$(basename "$i" | cut -d'-' -f1)
    gatk --java-options "-Xmx2g -Xms2g -XX:ParallelGCThreads=2" GenotypeGVCFs \
        -R ${reference} \
        -V gendb://${output_dir}/tmp/${interval_id} \
        -O ${output_dir}/${sample_name}_${interval_id}.vcf.gz
}
export -f run_GenotypeGVCFs
parallel --eta --jobs 3 run_GenotypeGVCFs ::: ${scattered_intervals[@]}
exit_value=$?
checkExitStatus "GenotypeGVCFs" $exit_value

#7a. Merge all 24 intervals into one vcf file, we use the Picard version of GatherVCFs, for simplicity this is called from gatk

gather_inputs="" #init variable to store list of scattered GVCFs
#Looping to append scattered GVCFs to gather_inputs if the file is found
for i in $(seq -w 0000 0023); do
    gather_inputs+="I=${output_dir}/${sample_name}_${i}.vcf.gz "
done

gather_inputs_len=$(echo $gather_inputs | wc -w) #as gather_inputs is a sequence of strings, use wc --words to count the number of entries
#error check to ensure files were correctly matched
if [ ${gather_inputs_len} -ne ${scatter_count} ]; then
    echo "ERROR: Mismatched of GVCF files, expected ${scatter_count}! Retry HaplotypeCaller."
    exit 1
else
    echo "File count OK: Found $gather_inputs_len GVCFs, running GatherVcfs..."
fi

gatk GatherVcfs \
    ${gather_inputs} \
    O=${output_dir}/${sample_name}.vcf.gz
exit_value=$?
checkExitStatus "GatherVcfs" $exit_value

# Remove intermediate files
echo "Removing intermediate files"
rm ${sortedMkdup_bam} ${sortedMkdup_bam}.bai ${sorted_bam} ${fastq1_qc} ${fastq2_qc}

for i in $(seq -w 0000 0023); do
    scattered_vcf=${output_dir}/${sample_name}_${i}.vcf.gz
    scattered_vcf_index=${output_dir}/${sample_name}_${i}.vcf.gz.tbi

    if [ -f ${scattered_vcf} ] && [ -f ${scattered_vcf_index} ]; then
        rm ${scattered_vcf} ${scattered_vcf_index}
    fi
done

checkScriptExit $exit_value
