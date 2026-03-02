#!/bin/bash
########################################################################
#.psmcfa-->.psmc with PSMC 
########################################################################
log_fold="path/to/logs"
psmc_fold="path/to/output/psmc_files"
mkdir -p "$log_fold" "$psmc_fold"

psmc="/path/to/psmc"

# List with .psmcfa files paths
file_list="$1"

while IFS= read -r file; do 
  ind=$(basename "$file" .psmcfa)
  ## Example with -p atomic interval = "2*3+58*1"
  ### Change -p according to the processed samples
  echo "[$(date)]  RUNNING $ind">>$log_fold/psmc_log.txt
  $psmc -N20 -t10 -r5 -p "2*3+58*1" -o $psmc_fold/$ind.psmc "$file"
  echo -e "[$(date)]\tPSMC\t$ind DONE">>$log_fold/psmc_log.txt
done<"$file_list"



