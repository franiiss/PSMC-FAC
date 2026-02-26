#!/bin/bash
##############################################################
#PREPARING .PSMCFA 
#bam-->vcf.gz+tabix-->fq.gz-->psmcfa
##############################################################
#Folder declaration
## Reference must be indexed with samtools faidx
ref="path/to/ref/genome.fa"
vcf_fold="/path/to/output/vcf"
fq_psmca_fold="/path/to/output/fq_psmca"
log_fold="path/to/logs"
mkdir -p "$vcf_fold" "$fq_psmca_fold" "$log_fold"
# Programs used
bcftools="path/to/bcftools"
htslib="/path/to/htslib"
progpsmc="/path/to/psmc"
vcfutils="/path/to/vcfutils.pl"

# List with all the paths to the original.bam files 
file_list="$1"

# Create dictionary with sample and coverage (from 0.depth_calc.sh)
declare -A dict_cvg
while IFS=, read -r sample coverage; do
    dict_cvg["$sample"]=$coverage
done < /path/to/samples_coverages.csv

while IFS= read -r file; do 
  # Cut the filename 
  ind=$(basename "$file" .bam)
  # Converting .bam into .vcf.gz and indexing
  echo  "[$(date)]-----START-----VCF-----  $ind">> $log_fold/psmcfa_log.txt
  $bcftools/bcftools mpileup -C50 -f $ref -Ou $file | $bcftools/bcftools call -c -Oz -o "$vcf_fold/$ind.vcf.gz"
  $htslib/tabix "$vcf_fold/$ind.vcf.gz"
  echo "[$(date)]-----DONE-----VCF----- $ind">> $log_fold/psmcfa_log.txt
  # Converting .vcf.gz --> .fq.gz  --> .psmcfa 
  ### Parameters -d is set to 5 and -D to x2 coverage mean
  cvg=${dict_cvg[$ind]}
  D_param=$(echo "$cvg * 2" | bc -l)
  
  echo "[$(date)]-----START-----PSMCFA----- $ind">> $log_fold/psmcfa_log.txt
  $bcftools/bcftools view "$vcf_fold/$ind.vcf.gz" | $vcfutils vcf2fq -d 5 -D $D_param | gzip > $fq_psmca_fold/$ind.fq.gz
  $progpsmc/utils/fq2psmcfa -q20 $fq_psmca_fold/$ind.fq.gz > $fq_psmca_fold/$ind.psmcfa
  echo "[$(date)]-----DONE-----PSMCFA----- $ind">> $log_fold/psmcfa_log.txt

# Delete previous files if not needed for other analysis
### Only leaves the .psmcfa file 
  rm $vcf_fold/$ind.vcf.gz
  rm $fq_psmca_fold/$ind.fq.gz
done < "$file_list"
  
