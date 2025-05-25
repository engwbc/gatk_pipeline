#!/usr/bin/env bash

SCRIPT_PATH=$(dirname "$(realpath "$0")") #Get script's path to source globalFunctions
source ${SCRIPT_PATH}/globalFunctions.sh #Import custom functions eg. createOutDir
exit_value= # Assign a variable for exit code $? 

output_as_bam=false #Initialise variable for only obtaining .bam file

while getopts "i:o:b:r:s:L:Bh" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;    #/path/to/fastq/folder
        o) output_dir="$OPTARG" ;;   #/path/to/output/folder
        b) bamFile="$OPTARG" ;;      #/path/to/input/bam
        r) reference="$OPTARG" ;;    #/path/to/reference/ref.fasta
        s) dbsnp="$OPTARG" ;;        #/path/to/knownSNPs.vcf.gz - for BaseRecalibration
        L) intervals="$OPTARG";;     #/path/to/variant_call_interval.file; in the future if scattered intervals exist, create a new option to skip SplitIntervals
        B) output_as_bam=true;;      #add -B as an argument to only generate a recalibrated bam file
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
if [[ -z ${bamFile} ]]; then
    checkInputArgs input_dir output_dir reference dbsnp intervals
else
    checkInputArgs bamFile output_dir reference dbsnp intervals
fi
log_file="fastq2vcf.log" #for record-keeping and debugging

#Check and create output directory if missing
createOutDir ${output_dir}

#Inititate logging
processLog ${output_dir}/${log_file} && start_time=$(logrunTime)
num_threads=$(nproc)
echo "[INFO] Number of CPU threads available: $num_threads"

#Check for dependencies
checkProgram bwa gatk samtools fastp java parallel

#For logging purposes
printf "Input directory: ${input_dir}\nReference: ${reference}\n\
Output directory: ${output_dir}\n"

# Initiate genomic region splitting for HaplotypeCaller
# this is done once per run as intervals can be reused
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

#Check if reference has been indexed, if not run bwa index
checkRefIndex $reference

# Define workflow functions
ApplyHaplotypeCaller() {
    haplo_outputfolder=${output_dir}/haplotypecaller
    createOutDir $haplo_outputfolder
    haplo_outputfile=${haplo_outputfolder}/${sample_name}

    echo "[INFO] Running HaplotypeCaller"
    # Export variables for parallel
    export bqsr_bam reference haplo_outputfile
    run_HaplotypeCaller() {
        i=$1
        interval_id=$(basename "$i" | cut -d'-' -f1)
        gatk --java-options "-Xmx8g" HaplotypeCaller \
            -mbq 20 \
            -I ${bqsr_bam} \
            -R ${reference} \
            -ERC GVCF \
            -L ${i} \
            -O ${haplo_outputfile}_${interval_id}.g.vcf.gz \
            --native-pair-hmm-threads 2
    }
    export -f run_HaplotypeCaller
    parallel --eta --jobs 3 run_HaplotypeCaller ::: "${scattered_intervals[@]}"
    checkExitStatus "HaplotypeCaller" $?
}

ApplyGenomicsDBImport() {
    ## Create a map file for batch importing into GenomicsDBImport to merge the g.vcf files
    echo "Found `find ${haplo_outputfolder} -type f -name "*\.g.vcf.gz" | wc -l` g.vcf.gz files" #find in haplotypecaller folder
    echo "Creating sample list for GenomicsDBImport"
    > ${output_dir}/${sample_name}.map #for debugging purposes the .map file is cleared if it exists from a previous run for whatever reason
    for file in ${haplo_outputfolder}/${sample_name}_*.g.vcf.gz; do
        interval_id=$(basename $file | sed -E 's/.*_([0-9]{4})\.g\.vcf\.gz/\1/') #extract interval id
        echo $interval_id
        map_sample_name=${sample_name}_${interval_id}
        echo -e "${map_sample_name}\t${file}" >> ${output_dir}/${sample_name}.map #write matched gvcf to map file
    done

    if [ -f "${output_dir}/${sample_name}.map" ]; then
        echo "[INFO] Finished creating ${sample_name}.map file"
    else
        echo "ERROR: Could not generate .map file for GenomicsDBImport!"
        exit 1
    fi    
    echo "[INFO] Running GenomicsDBImport"
    #need to create separate genDB folder for each interval... so ${output_dir}/genoDB/${interval}
    #based on NIH BioWulf benchmarks, increasing RAM does not improve performance
    #Sample GVCFs are compiled into a map file.
    export output_dir sample_name
    run_GenomicsDBImport() {
        i=$1
        interval_id=$(basename "$i" | cut -d'-' -f1)
        gatk --java-options "-Xmx2g -Xms2g -XX:ParallelGCThreads=2" GenomicsDBImport \
            --genomicsdb-workspace-path ${output_dir}/tmp/${sample_name}_${interval_id} \
            --sample-name-map ${output_dir}/${sample_name}.map \
            --tmp-dir ${output_dir}/tmp \
            --max-num-intervals-to-import-in-parallel 3 \
            -L ${i}
    }
    export -f run_GenomicsDBImport
    parallel --eta --jobs 3 run_GenomicsDBImport ::: ${scattered_intervals[@]}
    checkExitStatus "GenomicsDBImport" $?
}

ApplyGenotypeGVCFs() {
    echo "[INFO] Performing joint genotyping for $sample_name"
    echo "[INFO] For cohort studies (e.g. families), use CombineGCVFs prior to this step to compile multiple samples!"
    export reference output_dir sample_name
    run_GenotypeGVCFs() {
        i=$1
        interval_id=$(basename "$i" | cut -d'-' -f1)
        gatk --java-options "-Xmx2g -Xms2g -XX:ParallelGCThreads=2" GenotypeGVCFs \
            -R ${reference} \
            -V gendb://${output_dir}/tmp/${sample_name}_${interval_id} \
            -O ${output_dir}/${sample_name}_${interval_id}.vcf.gz
    }
    export -f run_GenotypeGVCFs
    parallel --eta --jobs 3 run_GenotypeGVCFs ::: ${scattered_intervals[@]}
    checkExitStatus "GenotypeGVCFs" $?
}

ApplygatherVCFs(){
    # Merge all 24 intervals into one vcf file using the Picard version of GatherVCFs, for simplicity this is called from gatk
    #init variable to store list of scattered GVCFs
    gather_inputs=""
    #Looping to append scattered GVCFs to gather_inputs if the file is found 
    for i in $(seq -w 0000 0023); do
        gather_inputs+="I=${output_dir}/${sample_name}_${i}.vcf.gz "
    done

    #as gather_inputs is a sequence of strings, use wc --words to count the number of entries
    gather_inputs_len=$(echo $gather_inputs | wc -w)
    if [ ${gather_inputs_len} -ne ${SCATTER_COUNT} ]; then
        echo "ERROR: Mismatched number of GVCF files, expected ${SCATTER_COUNT}! Retry HaplotypeCaller."
        exit 1
    else
        echo "File count OK: Found $gather_inputs_len GVCFs, running GatherVcfs..."
    fi

    gatk GatherVcfs \
        ${gather_inputs} \
        O=${output_dir}/${sample_name}.vcf.gz
    exit_value=$?
    checkExitStatus "GatherVcfs" $exit_value
    #-I input.001.vcf.gz -I input.002.vcg.gz -O output.vcf.gz
    # include java -Dpicard.useLegacyParser=false
}

removeTempFiles() {
    local tempFiles=()
    for args in "$@"; do
        tempFiles+=("$args")
    done
    for file in "${tempFiles[@]}"; do
        if [ -f "$file" ]; then
            rm -f "$file"
            echo "[INFO] Removed: $file"
        else
            echo "[WARNING] Could not find: $file"
        fi
    done
}

# Separate pipeline into two branches: a) normal run (fastq to vcf) and b) bam to vcf
# Branched pipeline will merge at the Variant Calling step
workflowMain() {
    if [[ -n "${input_dir}" && -z "${bamFile}" ]]; then
        while read fastq; do
            ## Retrieve fastQ files ##
            echo "[INFO] Searching for samples in ${input_dir}."
            echo "[INFO] find command is set to -maxdepth 1. To avoid errors, provide the full path to the FOLDER containing the fastQ file(s)."
            sample_StartRunTime=$(logrunTime)
            sample_name=$(basename "$fastq" | sed 's/_1.*$//')
            fastq1=$fastq #Takes read 1 as either fq.gz/fastq.gz
            #If file extension naming is inconsistent, the following code block checks and allows only one extension type.
            fastq2_1=${input_dir}/${sample_name}_2.fq.gz
            fastq2_2=${input_dir}/${sample_name}_2.fastq.gz

            if [ -f "$fastq1" ]; then
                echo "Found read 1: $fastq1"
            else
                echo "ERROR: Read 1 for ${sample_name} not found!"
                end_time=$(logrunTime) && endLog "$start_time" "$end_time"
                exit 2
            fi

            fastq2_list=()
            fastq2_list=($(checkRead2File "$fastq2_1" "$fastq2_2"))
            #From fastq2_list, index the valid file from the array 
            fastq2="${fastq2_list[0]}"
            if [ -f "$fastq2" ]; then #second sanity check 
                echo "Found read 2: $fastq2"
            else
                echo "ERROR: read 2 for ${sample_name} not found!"
                end_time=$(logrunTime) && endLog "$start_time" "$end_time"
                exit 2
            fi

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
            #CPU usage: peaks at ~40% on 4 threads (Xeon @2.1Ghz)
            exit_value=$?
            checkExitStatus "fastp" $exit_value
            echo "[INFO] Please check fastp output with multiqc!"

            fastq1_qc=${QC_DIR}/${sample_name}_qc_1.fq.gz
            fastq2_qc=${QC_DIR}/${sample_name}_qc_2.fq.gz

            ################################################
            ### Step 2. Alignment and positional sorting ###
            ################################################

            echo "Starting BWA-MEM alignment and sorting for ${sample_name}"
            printf "\n[INFO] Normal bwa mem stdout is [M::mem_pestat]. \
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

            if [ "$output_as_bam" = true ]; then
                rm -r ${output_dir}/tmp
                rm ${sortedMkdup_bam} ${sortedMkdup_bam}.bai ${sorted_bam} ${fastq1_qc} ${fastq2_qc}
                checkScriptExit $exit_value
                exit 0
            fi
            
            # GATK HaplotypeCaller #

            ApplyHaplotypeCaller

            # GenotypeGVCFs (gvcf to vcf) #

            ApplyGenomicsDBImport
            ApplyGenotypeGVCFs

            #  Merge gVCFs    #

            ApplygatherVCFs

            # Remove intermediate files
            removeTempFiles ${sortedMkdup_bam} ${sortedMkdup_bam}.bai ${sorted_bam} ${fastq1_qc} ${fastq2_qc}
            # Remove scattered vcfs and their indices 
            for i in $(seq -w 0000 0023); do
                scattered_vcf=${output_dir}/${sample_name}_${i}.vcf.gz
                scattered_vcf_index=${output_dir}/${sample_name}_${i}.vcf.gz.tbi
                removeTempFiles ${scattered_vcf} ${scattered_vcf_index}
            done
            
        done < <(find ${input_dir} -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz" \))
    elif [[ -n "${bamFile}" ]]; then
        # If option -b is given: this code block is run.
        sample_StartRunTime=$(logrunTime)
        echo "[INFO] Option -b given, converting bam to vcf"
        sample_name=$(basename ${bamFile} | cut -d'_' -f1-3)
        echo "[INFO] Input file: ${bamFile}"
        echo "[INFO] Sample name: ${sample_name}"
        bqsr_bam=${bamFile} # this assumes the input bam file has been base recalibrated!
        samtools index ${bqsr_bam} # re-indexes the bam file just to be safe
        exit_value=$?
        checkExitStatus "Index BQSR" $exit_value
        
        ################################
        # Step 6. GATK HaplotypeCaller #
        ################################

        ApplyHaplotypeCaller

        ###########################################
        #   Step 7. GenotypeGVCFs (gvcf to vcf)   #
        ###########################################

        ApplyGenomicsDBImport
        ApplyGenotypeGVCFs

        #############################
        #   Step 7.5 Merge gVCFs    #
        #############################

        ApplygatherVCFs

        # Remove scattered vcfs and index files
        for i in $(seq -w 0000 0023); do
            scattered_vcf=${output_dir}/${sample_name}_${i}.vcf.gz
            scattered_vcf_index=${output_dir}/${sample_name}_${i}.vcf.gz.tbi
            removeTempFiles ${scattered_vcf} ${scattered_vcf_index}
        done
    else
        echo "[ERROR] No valid input provided!"
        exit 1
    fi
}

# Run pipeline
workflowMain
rm -r ${haplo_outputfolder}
rm -r ${output_dir}/tmp

#End logging
sample_endRunTime=$(logrunTime) && endLog ${sample_StartRunTime} ${sample_endRunTime} ${sample_name}
checkScriptExit $exit_value
