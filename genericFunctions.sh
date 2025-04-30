#!/usr/bin/env bash

printUsage() {
    printf "\nUsage: bash /path/to/${0} -i /path/to/fastq/folder -o /path/to/output/folder -r /path/to/reference.fasta -s /path/to/snpDB.vcf -L /path/to/your.intervals_list\n"
}

#Check if required programs are installed
checkProgram() {
    local missing_programs=() # Array to store missing programs
    for program in "$@"; do
        #use bash command with argument -v to execute the programs by searching system PATH
        #redirect stdout and stderror to /dev/null to proceed to the conditional check 
        if ! command -v "$program" >/dev/null 2>&1; then   
            missing_programs+=("$program") # Add missing program to the array
        fi
    done
    if [ ${#missing_programs[@]} -ne 0 ]; then
        echo "ERROR: Could not find the following tools in your environment:"
        for program in "${missing_programs[@]}"; do
            echo "  - $program"
        done  
        echo "Please try installing them using your package manager (e.g., apt, yum, conda, etc.)"
        exit 1
    fi
}

#Check input directory and reference files
checkInputArgs(){
    local required_options=("input_dir" "output_dir" "reference" "dbsnp" "intervals")
    
    for option in ${required_options[@]}; do
        if [[ -z ${!option} ]]; then
            echo "ERROR: Missing required arguments!"
            printUsage "$@"
            exit 1
        fi
    done
}

createOutDir(){
    local dir=$1
    if [ ! -d $dir ]; then
      echo "Creating output folder: $dir"
      mkdir -p $dir #create parent directories if needed, but should be somewhere the user has permission. 
      echo "Created $dir"
    else
        echo "Output directory already exists"
    fi
}

createTmpDir(){
    local tmp_dir=$1
    if [ ! -d $tmp_dir ]; then
      echo "Creating temp folder: $tmp_dir"
      mkdir -p $tmp_dir
    else
      echo "$tmp_dir exists"
    fi
}

#Save stdout to a log file
processLog(){
   
    log_file=$1
    base_filename="${log_file%.*}"
    file_ext="${log_file##*.}" #Get file extension (e.g., .txt or .log)
    
    counter=0
    while [ -f ${log_file} ]; do
      if [ $counter -eq 0 ]; then
        log_file=${base_filename}.${file_ext}
      else
        log_file=${base_filename}${counter}.${file_ext}
      fi
      counter=$((counter+1))
    done
    
    touch "$log_file" #Create log file if it doesn't exist
    echo "Logging to: $log_file"
    start_time=$(date +%s)
    exec > >(tee -a "$log_file") 2>&1
    echo "Log time: `date -d @$(date +%s)`" #Echo human-readable date and time
    echo "Using bash from: `which bash`" #Check if right interpreter is called
    #Checks if script was cancelled during logging"
    trap 'end_time=$(date +%s); runtime=$((end_time - start_time)); echo "Script was terminated by user at $(date -d @$end_time)."; echo "Runtime before termination: $runtime seconds."; exit 1' SIGINT
}

#Get time for logging
logrunTime(){
  time=$(date +%s)
  echo "$time"
}

#Calculate total script runtime
endLog(){
  local start_time=$1
  local end_time=$2
  local sample=$3
  runtime=$((end_time - start_time))
  if [ -z ${sample} ]; then
    echo "Total runtime: ${runtime} seconds"
  else
    echo "${sample} runtime: ${runtime} seconds"
  fi
  #exec > /dev/tty #revert stdout to default terminal
}

#Exit status
#checkExitStatus "Alignment" $exit_value
checkExitStatus(){
  local process_name=$1
  local exit_value=$2
  if [ $exit_value -ne 0 ]; then
      printf "ERROR %d: Could not perform %s\n" "$exit_value" "$process_name"
      end_time=$(logrunTime) && endLog "$start_time" "$end_time"
      exit $exit_value
  else
      printf "%s: OK\n" "$process_name"
  fi
}

#Function to check whether read2.fq.gz or read2.fastq.gz exists for read 2; will only accept one - whichever exists but not both for naming consistency and file compatibility.
checkRead2File(){
  local fastq2_1="$1"
  local fastq2_2="$2"
  #Append possible filenames to array
  fastq2_list+=("$fastq2_1")
  fastq2_list+=("$fastq2_2")
  fastq2_list=($(for file in "${fastq2_list[@]}"; do #iterate through array and if the fastq file exists, we keep it.
    if [[ -f "$file" ]]; then
      echo "$file"
    fi
  done))
  #To make troubleshooting easy, we notify the user that the file is missing or multiple files are found.
  if [ ${#fastq2_list[@]} -eq 0 ]; then
      echo "ERROR: No valid FASTQ files for read 2 found."
      end_time=$(logrunTime) && endLog "$start_time" "$end_time"
      exit 2
  elif [ ${#fastq2_list[@]} -gt 1 ]; then
      echo "ERROR: Expected one FASTQ file. Found ${#fastq2_list[@]} files."
      end_time=$(logrunTime) && endLog "$start_time" "$end_time"
      exit 2
  fi
  
  echo "${fastq2_list[@]}" #return list to fastq2 variable in the main script
}

#Check if the reference genome has already been indexed
checkRefIndex(){
  local reference=$1
  local ref_name=$(basename ${reference})
  local ref_path=$(dirname ${reference})
  #bwa index produces 5 output files; these are listed in an array, index_files
  local index_files=("${ref_path}/${ref_name}.amb" "${ref_path}/${ref_name}.ann" \
  "${ref_path}/${ref_name}.pac" "${ref_path}/${ref_name}.sa" "${ref_path}/${ref_name}.bwt")
  local missing_index=()
  #Add missing index files to array
  for file in ${index_files[@]}; do
    if [ ! -f $file ]; then
      missing_index+=($file)
    fi
  done

  #If any index files are missing, re-index the ref.
  if [ ${#missing_index[@]} -gt 0 ]; then
    echo "Detected missing bwa files: ${missing_index[*]}"
    echo "Running bwa index for: $ref_name"
    bwa index ${reference}
    exit_value=$?
    checkExitStatus "Index reference" $exit_value
  else
    echo "${ref_name} has been indexed, proceeding..."
  fi
}

#Final check to ensure smooth exit + end logging
checkScriptExit(){
  local exit_value=$1
  if [ $exit_value -eq 0 ]; then
    echo "Exit: $exit_value"
    echo "Completed: `date -d @$(date +%s)`"
    end_time=$(logrunTime) && endLog "$start_time" "$end_time"
  else
    echo "Pipeline was terminated due to an error. Exit: $exit_value"
      end_time=$(logrunTime) && endLog "$start_time" "$end_time"
      exit $exit_value
  fi 
}