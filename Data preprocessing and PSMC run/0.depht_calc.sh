#!/bin/bash

out_folder="/path/to/output/folder"
log_folder="/path/to/log/folder"
# A text file ,which contains the paths of the files to be analyzed (1 per line),is inserted as the first argument
## Example: bash 0.depth_calc.sh /path/to/desired/list_file
file_list="$1"

# Create a folder if needed
mkdir -p $out_folder/
# Activate conda environment with samtools (if needed) 
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate samtools

# Iterate through the lines of the file_list 
## Call samtools depth for every path/line in that list
### Create 2 files: depths_file will store the coverage and standard deviation 
### 		    csv file will be used to create a dictionary in the downsampling step
while IFS= read -r input_file; do 
	ind=$(basename "$input_file" .autosomes.bam)
	echo -e "[$(date)]\tSTARTED-------COVERAGE-ANALYSIS-----$ind">>$log_folder/depth_log.txt
	samtools depth -a "$input_file" | awk -v bam="$input_file" -v depths="$out_folder/depths_file" -v csv="$out_folder/coverages.csv" '{sum+=$3; sumsq+=$3*$3} END {print bam  "      Average : " sum/NR "      Stdev = " sqrt((sumsq-sum^2/NR)/NR) >> depths; print bam","sum/NR >> csv}'
	echo -e "[$(date)]\tFINISHED-------COVERAGE-ANALYSIS------$ind">>$log_folder/depth_log.txt
done < "$file_list"

conda deactivate


