# PSMC-FAC
False-negative rate Automated Correction in low-coverage genomes: A statistical framework for correcting loss of heterozygosity in low-coverage genomic demographic inference


This repository compiles the scripts needed to perform PSMC analysis coupled with FNR correction.

Contents:

- [Initial Calculations](#initial-calculations)
    - [Depth](#depth)
    - [Heterozygosity](#heterozygosity)
- [Downsampling](#downsampling)
- [File conversion](#file-conversion)
- [PSMC run](#psmc-run)

- [Distance calculation](#distance-calculation)

- [Regression analysis](#regression-analysis) 
- [Plot samples](#plot-samples)
- [Tips And Considerations](#tips-and-considerations)


# Initial Calculations
To asses the charachteristics of the samples used, check the coverage (read depth) and heterozygosity. 
## Depth
samtools depth -a (all positions) piped to awk was used to output the average coverage and the standard deviation.
```bash
samtools depth -a "$input_file" | awk -v bam="$input_file" '{sub(/.*\//, "", bam); sum+=$3; sumsq+=$3*$3} END {print bam  "      Average : " ,sum/NR "      Stdev = ", sqrt((sumsq-sum^2/NR)/NR)}' >> $folder/depthsfile
```
## Heterozygosity
Making use of ANGSD and samtools the bash script does the following:

 FOR WOLF FOR NOW 38chr... (lo cambiaré)
- Extracts every chromosome in the bam file
```bash
for i in {1..38}; do
        chr=$(printf "%02d" "$i")
        samtools view -b $file chr$chr > $new_folder/$ind.chr$chr.bam
	    samtools index $new_folder/$ind.chr$chr.bam
done
```

- Creates a .sfs file per .bam file extracted. 
	- First, ANGSD is called to process the bam file.
 	- Then, realSFS function from ANGSD processes the .saf.idx file outputed by ANGSD to generate the unfolded .sfs file.
```bash
for i in {1..38}; do
        chr=$(printf "%02d" "$i")
        $angsd -i "$new_folder/$ind.chr$chr.bam" -doSaf 1 -out $new_folder/angsd/chr$chr -anc $ref/canfam31.chr$chr.fa -GL 2 -P 4 -minQ 20 -minMapQ 20 
        $real $new_folder/angsd/chr$chr.saf.idx -P 4  > $new_folder/angsd/SFS/chr$chr.sfs
done
```
Lastly, the unfolded .sfs file contains the iformation of  Homozygous (0/0), Heterozygous (0/1) and Homozygous (1/1) sites in this order.
The Heterozygosity equals dividing the Heterozygous by the total sites. 
```bash
col1=0
col2=0
col3=0
for i in {1..38}; do
    chr=$(printf "%02d" "$i")
    sfs_file="$new_folder/angsd/SFS/chr$chr.sfs" 
    # Read and add values from each SFS file
    while read -r c1 c2 c3; do
        col1=$(echo "$col1 + $c1" | bc)
        col2=$(echo "$col2 + $c2" | bc)
        col3=$(echo "$col3 + $c3" | bc)
    done < "$sfs_file"
done
heterozygosity1=$(awk '{hets=$2; totsites=$1+$2+$3; heterozygosity=hets/totsites; print heterozygosity}' "$new_folder/angsd/SFS/combined.$ind.sfs")
```
# Downsampling
The first step of the FNR calculation is downsampling a high-coverage sample to a determine set of coverage values.<br/>
> samtools view -h (include header) -s [int value (0,1)]

s determines the fraction of reads kept in the output file
```bash
samtools view -h -s 0.6756 -b $file > $downfold/$ind.15.bam
```
In this example, by retaining 67.56% reads of the original bam file (22.202X) the output bam file will present ~15X coverage.

# File conversion
To prepare the .psmcfa files needed for PSMC, variant calling on the bam files is performed with bcftools mpileup and call functions. Indexing the .vcf.gz output file at the end.
```bash
$bcftools/bcftools mpileup -C50 -f $ref -Ou $file | $bcftools/bcftools call -c -Oz -o "$vcfold/$filename.vcf.gz"
$htslib/tabix "$vcfold/$filename.vcf.gz"
```
Next, convert the .vcf files to .fq format using bcftools view of the .vcf.gz file piped to vcf2fq -d -D piped to gzip to obtain the .fq.gz file.

During this step, the next vcf2fq parameters need to be determined:

> -d sets the lower coverage threshold, set to 5.
> D sets the higher coverage threshold, set to 2xAvgCvg (following recommendation from the official PSMC github)
```bash
# coverage of every file is stored in an associative array
cvg=$(dict_cvg[$filename])
D_param=$(echo "$cvg * 2" | bc -l)
$bcftools/bcftools view "$vcfold/$filename.vcf.gz" | $bcftools/bin/vcfutils.pl vcf2fq -d 5 -D $D_param | gzip > $fq_psmcafold/$filename.fq.gz
```
Lastly call fq2psmcfa (from /psmc/utils/ folder) to convert the .fq.gz file to .psmcfa.
```bash
$progpsmc/utils/fq2psmcfa -q20 $fq_psmcafold/$filename.fq.gz > $fq_psmcafold/$filename.psmcfa
```

# PSMC run
Pairwise Sequentially Markovian Coalescent (PSMC; Li & Durbin,2011) is called on the .psmcfa files obtained previously to output the .psmc file which contains the demographic information of the sample. To include the effect of FNR correction during PSMC's course two options are available.
 
PSMC can run on the selected file/s with -r fixed to 5. PSMC run example with four different -p settings. 
```
$psmc -N20 -t10 -r5 -p "1*6+58*1" -o $plotf/1_6/$filename.psmc "$file"
$psmc -N20 -t10 -r5 -p "2*3+58*1" -o $plotf/2_3/$filename.psmc "$file"
$psmc -N20 -t10 -r5 -p "3*2+58*1" -o $plotf/3_2/$filename.psmc "$file"
$psmc -N20 -t10 -r5 -p "3+3+58*1" -o $plotf/3__3/$filename.psmc "$file"
```
> [!NOTE]
> You can erase the PSMC command with -p "1 * 6 + 58 * 1" as it was only kept to prove Nadachowska's work.
> The results from the different -p used is quite similar, to reduce the processing time just run with one command.
> Also, take into account -p setting depends on multiple elements. You can make use of proven atomic intervals in trusted bibliography for the species studied.


## Regression analysis
Repeating the 2 plotting steps with all the downsamples generates a table linking downsample coverage and FNR value required for correction. 

<ins>A polynomial regression of degree 2</ins> can be fitted to this data, showing a high coefficient of determination (this method is still under examination and refinement).


# Plot samples
Once the quadratic equation from the regression analysis was accomplished, the samples below 18X were assigned a FNR value dependant of their coverage. This threshold was applied following Nadachowska et al. (2016) recommendation for samples with enough coverage (>18X) for PSMC analysis.

## References
Conda was used to create controlled environments to run specific tools. Conda documentation. (2024). https://docs.conda.io/


Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O.,
Whitwham, A., Keane, T., Mccarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years
of SAMtools and BCFtools. 10, 1–4. https://doi.org/10.1093/gigascience/giab008

Korneliussen, T. S., Albrechtsen, A., & Nielsen, R. (2014). ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics, 15(1), 1–13.
https://doi.org/10.1186/S12859-014-0356-4/TABLES/4

Li, H., & Durbin, R. (2011). Inference of human population history from individual
whole-genome sequences. Nature 2011 475:7357, 475(7357), 493–496.
https://doi.org/10.1038/nature10231


Hilgers, L., Liu, S., Jensen, A., Brown, T., Cousins, T., Schweiger, R., Guschanski,
K., & Hiller, M. (2025). Avoidable false PSMC population size peaks occur across
numerous studies. Current Biology, 35(4), 927-930.e3.
https://doi.org/10.1016/J.CUB.2024.09.028

Nadachowska-Brzyska K, Burri R, Smeds L, Ellegren H. PSMC analysis of effective population sizes in molecular ecology and its application to black-and-white Ficedula flycatchers. Mol Ecol. 2016 Mar;25(5):1058-72. doi: 10.1111/mec.13540. Epub 2016 Feb 15. PMID: 26797914; PMCID: PMC4793928.







+Contact information: fran.iglesiasantos@gmail.com
