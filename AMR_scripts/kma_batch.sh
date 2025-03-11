#!/bin/bash

# Start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Batch Started at: $start_time"

# specify directories (modify these paths before running)
in_dir="/path/to/raw_data"
out_dir="/path/to/output"
main_dir="/path/to/main_dir"
index_dir="/path/to/kma_index/"
batch_file="${main_dir}/samples.txt"

# create or clean the log file
> "${main_dir}/kma_log.txt"

# initialize tracker
sample_counter=0

# loop through all samples
while read -r line; do

   # Check if this is read1 or read2
   if [[ $line =~ \.1\.fastq\.gz$ ]]; then
    	
      ((sample_counter++))  # Increment sample counter
    	
      # Start time
      start_time2=$(date +"%Y-%m-%d %H:%M:%S")
      echo " ~~~ Sample "$sample_counter" Started at: $start_time2 ~~~ "

      echo " ** Extracting prefix and samples file name** "
    	
      # Extract prefix from the line
      prefix=$(echo "$line" | cut -d'_' -f2)
    	
      # It's read1, construct the command with read1 and read2
      read1=$(grep "${prefix}.*\.1\.fastq\.gz$" "$batch_file")
      read2=$(grep "${prefix}.*\.2\.fastq\.gz$" "$batch_file")
        	
      # Print a message for logging or verification
      echo ">> Sample: $prefix"
      echo ">> Read 1 = $read1"
      echo ">> Read 2 = $read2"

      # Run KMA
      ./kma -ipe "${in_dir}/${read1}" "${in_dir}/${read2}" -o "${out_dir}/${prefix}" -t_db "$index_dir" -1t1 -nc -na -nf -ef -ID 80 -ml 60
        
      # End time
      end_time2=$(date +"%Y-%m-%d %H:%M:%S")
      echo -e " ~~~ Sample "$sample_counter" Ended at: $end_time2 ~~~ \n"
   fi

done < "$batch_file" >> "${main_dir}/kma_log.txt" 2>&1

# End time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Batch Ended at: $end_time"
