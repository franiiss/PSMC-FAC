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
mkdir -p $vcf_fold $fq_psmca_fold $log_fold
# Programs used
bcftools="path/to/bcftools"
htslib="/path/to/htslib"
progpsmc="/path/to/psmc"
vcfutils="/path/to/vcfutils.pl"
# List with all the bams downsamples paths
file_list="$1"
while IFS= read -r file; do 
  ind=$(basename "$file" .bam)
  # Converting .bam into .vcf.gz and indexing
  echo "[$(date)]    Starting VCF of $ind">> $log_fold/psmcfa_log.txt
  $bcftools/bcftools mpileup -C50 -f $ref -Ou $file | $bcftools/bcftools call -c -Oz -o "$vcf_fold/$ind.vcf.gz"
  $htslib/tabix "$vcf_fold/$ind.vcf.gz"
  echo "[$(date)]    VCF of $ind done">> $log_fold/psmcfa_log.txt
  
  # Converting .vcf.gz --> .fq.gz  --> .psmcfa 
  ### Parameters -d is set to 5 and -D to x2 coverage mean
  ### Regular expression to obtain the coverage from the filename, will pick the number at the end just before ".bam"
  #### Filename must have a number indicating its coverage at the end
  #### Example: human243.9.bam --> X9
  cvg=$(echo $ind | grep -oE '[0-9]+$')
  D_param=$(( cvg * 2 ))

  echo "[$(date)]    Starting PSMCFA-ING OF  $ind">> $log_fold/psmcfa_log.txt
  $bcftools/bcftools view "$vcf_fold/$ind.vcf.gz" | $vcfutils vcf2fq -d 5 -D $D_param | gzip > $fq_psmca_fold/$ind.fq.gz
  $progpsmc/utils/fq2psmcfa -q20 $fq_psmca_fold/$ind.fq.gz > $fq_psmca_fold/$ind.psmcfa
  echo "[$(date)]    PSMCFA-ING OF  $ind done">> $log_fold/psmcfa_log.txt

  # Delete previous files if not needed
    # Left with the .psmcfa file
  rm $vcf_fold/$ind.vcf.gz
  rm $fq_psmca_fold/$ind.fq.gz
done < "$file_list"
