#!/bin/bash
#################################
# DOWNSAMPLE ALL FILES IN THE FILE-LIST
#################################
out_fold="/path/to/downsample_bams"
log_fold="/path/to/logs"
file="$1"

# After having the cvg calculated for all samples create a dictionary with sample_name & coverage
### With the csv file created with 0.depth_calc.sh
declare -A dict_cvg
while IFS=, read -r sample coverage; do
    dict_cvg["$sample"]=$coverage
done < /path/to/coverages.csv

# Activate conda enviroment with samtools (if needed)
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate samtools

# Adjust the last element to your own filename to leave just the sample name
ind=$(basename $file .bam)
orig_cvg=${dict_cvg["$ind"]}
  
echo -e "[$(date)]\t DOWNSAMPLING \t$ind" >>$log_fold/downsample_log.txt
# Calculate the downsample factor for the final wanted coverage (15,14,13,12,...,2)
for target_cvg in 15 14 13 12 11 10 9 8 7 6 5 4; do
	# The downsampling factor= -s flag, sets the % (-s 0,5 -> 50%) of reads kept in the output
	echo "[$(date)]\t Started Downsampling  $ind  to Coverage= $target_cvg" >>$log_fold/downsample_log.txt
	ratio=$(echo "scale=6; $target_cvg / $orig_cvg" | bc)
	output_file="$out_fold/$ind.$target_cvg.bam"
	samtools view -h -s "$ratio" -b "$file" > "$output_file"
	echo -e "[$(date)]\t$ind Downsample to Coverage= $target_cvg \t COMPLETED" >>$log_fold/downsample_log.txt
done
echo -e "[$(date)]\tFINISHED DOWNSAMPLING OF\t$ind" >>$log_fold/downsample_log.txt

conda deactivate
