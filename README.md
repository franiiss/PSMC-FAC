
This repository compiles the scripts needed to perform PSMC analysis coupled with FNR correction.

Contents:
- [Preprocessing](#preprocessing)
	- [Depth](#depth)
	- [Downsampling](#downsampling)
	- [File conversion](#file-conversion)
	- [PSMC run](#psmc-run)
- [PSMC-FAC](#psmc-fac)

# Preprocessing 


## Depth
samtools depth -a (all positions) piped to awk was used to output the average coverage and the standard deviation in a file.
A .csv file is also created to store the filename and coverage to create a dictionary during the downsampling step

```bash
samtools depth -a "$input_file" | awk -v bam="$input_file" -v depths="$out_folder/depths_file" -v csv="$out_folder/coverages.csv" '{sum+=$3; sumsq+=$3*$3} END {print bam  "      Average : " sum/NR "      Stdev = " sqrt((sumsq-sum^2/NR)/NR) >> depths; print bam","sum/NR >> csv}'
```
## Downsampling
The first step of the FNR calculation is downsampling a high-coverage sample to a determine set of coverage values.<br/>
> samtools view -h (include header) -s [int value (0,1)]

s determines the fraction of reads kept in the output file
```bash
samtools view -h -s 0.6756 -b $file > $downfold/$ind.bam
```
In this example, by retaining 67.56% reads of the original bam file (22.202X) the output bam file will present ~15X coverage.<br/>

The script interates through a set of desired coverages and calculates the -s ratio for each coverage dividing the target coverage by the original coverage.
```bash
for target_cvg in 15 14 13 12 11 10 9 8 7 6 5 4; do
	ratio=$(echo "scale=6; $target_cvg / $orig_cvg" | bc)
	output_file="$out_fold/$ind.$target_cvg.bam"
	samtools view -h -s "$ratio" -b "$file" > "$output_file"
done
```

## File conversion
To prepare the .psmcfa files needed for PSMC, variant calling on the bam files is performed with bcftools mpileup and call functions. Indexing the .vcf.gz output file at the end.
```bash
  $bcftools/bcftools mpileup -C50 -f $ref -Ou $file | $bcftools/bcftools call -c -Oz -o "$vcf_fold/$ind.vcf.gz"
  $htslib/tabix "$vcf_fold/$ind.vcf.gz"
```
Next, convert the .vcf files to .fq format using bcftools view on the .vcf.gz file piped to vcf2fq piped to gzip to obtain the .fq.gz file.

During this step, the next vcf2fq parameters need to be determined:

> -d sets the lower coverage threshold, set to 5 <br/>
> -D sets the higher coverage threshold, set to 2xAvgCvg (following recommendation from the official PSMC github)

<br/>
How the calculation of -D values is performed varies between the original bam files and the downsamples. The coverage of original bam files was stored in a .csv with 0.depth_calc.sh and it's loaded in a dictionary, while for the downsamples a simple regression expression extracts the coverage from the filename.


```bash
# For Downsamples
  cvg=$(echo $ind | grep -oE '[0-9]+$')
# For Original files
  cvg=${dict_cvg[$ind]}

### Both share the next steps
  D_param=$(( cvg * 2 ))
  $bcftools/bcftools view "$vcf_fold/$ind.vcf.gz" | $vcfutils vcf2fq -d 5 -D $D_param | gzip > $fq_psmca_fold/$ind.fq.gz
```
Lastly call fq2psmcfa (from /psmc/utils/ folder) to convert the .fq.gz file to .psmcfa.
```bash
  $progpsmc/utils/fq2psmcfa -q20 $fq_psmca_fold/$ind.fq.gz > $fq_psmca_fold/$ind.psmcfa
```

## PSMC run
Pairwise Sequentially Markovian Coalescent (PSMC; Li & Durbin,2011) is called on the .psmcfa files obtained previously to output the .psmc files.<br/>
>-p: time_vector (species-specific)<br/>
>-N: maximum number of iterations<br/>
>-t: maximum 2N0 coalescent time<br/>
>-r: initial theta/rho ratio
```
  $psmc -N20 -t10 -r5 -p "2*3+58*1" -o $psmc_fold/$ind.psmc "$file"
```

> [!NOTE]
> You can erase the PSMC command with -p "1 * 6 + 58 * 1" as it was only kept to prove Nadachowska's work.




# PSMC-FAC

**PSMC-FAC** False-negative rate Automated Correction in low-coverage genomes: A statistical framework for correcting loss of heterozygosity in low-coverage genomic demographic inference. An unified Python workflow designed to evaluate and optimize false-negative-rate (FNR) corrections for Pairwise Sequentially Markovian Coalescent (PSMC) trajectories inferred from high- and reduced-coverage genomes.

The pipeline integrates:

1. Extraction of FNR-corrected demographic trajectories  
2. Computation of Hausdorff and discrete Fréchet distances  
3. Calculation of residual errors (SSE and log10-SSE)  
4. Identification of optimal FNR values across coverage levels  
5. Optional polynomial fitting of coverage–FNR relationships  
6. Optional visualization of corrected demographic trajectories  

All analyses are executed from a single command-line interface.

---

# Summary

For each downsampled genome, PSMC-FAC evaluates a discrete grid of candidate FNR values spanning:

```
[0, 0.99]
```

in increments of 0.01 (configurable).

Each candidate FNR produces a corrected PSMC trajectory.  
Thus, each downsampled genome generates up to **100 corrected demographic curves**.

These corrected curves are compared to the corresponding high-coverage reference trajectory.

---

# Mathematical Formulation

All trajectories are projected onto a shared logarithmically spaced temporal grid.

Each curve is represented as:


$$P   = {(t_i, Ne_i^ref)}_{i=1}^n$$
$$Q_f = {(t_i, Ne_i^(f))}_{i=1}^n$$


Where:

- `P`   = reference trajectory  
- `Q_f` = FNR-corrected trajectory  
- `t_i` = shared time grid  
- `Ne_i` = effective population size  

Distances computed:

- Symmetric Hausdorff distance  
- Discrete Fréchet distance  

Residual errors computed:

- SSE (linear scale)  
- SSE (log10 scale)

---

# Workflow Overview

PSMC-FAC performs five integrated steps:

## 1. Extract FNR-corrected trajectories
- Reads final iteration of each `.psmc`
- Applies FNR grid
- Resamples on custom log-spaced time vector
- Optionally writes one TSV per PSMC file

## 2. Compute geometric distances
- Hausdorff
- Fréchet

Distances computed under user-defined coordinate scaling:

- `psmc` (default; log10 time + Ne × 10⁴)
- `linear`
- `loglog`

## 3. Compute residual errors
- SSE
- SSE_log10

## 4. Identify optimal FNR
For each sample × coverage:
- FNR minimizing Hausdorff
- FNR minimizing Fréchet

## 5. Optional polynomial fitting
Fits:

```
FNR = ax² + bx + c
```

Where:
- `x` = coverage
- `y` = optimal FNR

(Only performed if ≥ 3 coverage levels)

## 6. Optional trajectory plots
Generates multi-panel PDFs comparing:

- Reference
- Uncorrected (FNR=0)
- Hausdorff-optimal
- Fréchet-optimal

---

# Directory Structure Assumption

Within each final directory level:

```
sample_directory/
├── sample.psmc          (reference)
├── sample.4.psmc        (4x coverage)
├── sample.5.psmc
├── sample.6.psmc
```

Reference and downsampled files **must reside in the same directory**.

---

# Example Usage (Wolf Example from Repository)

This is the exact example provided in:

`example_run_psmc_fac_wolf.sh`

```bash
python3 PSMC-FAC.py \
  --root example_run_PSMC_FAC_wolf \
  --base_files gr_wolf.spain.SAMN04851099.psmc \
  --outdir example_run_PSMC_FAC_output_wolf \
  --mu 4.5e-9 \
  --g 4.4 \
  --FNR_min 0.0 \
  --FNR_max 0.99 \
  --svalue 100 \
  --tmin 1e4 \
  --tmax 1.5e6 \
  --write_psmc_tsv
```

This command:

- Recursively scans `example_run_PSMC_FAC_wolf`
- Treats `gr_wolf.spain.SAMN04851099.psmc` as the baseline
- Evaluates all other `.psmc` files as downsampled genomes
- Tests FNR values from 0.0 to 0.99
- Writes corrected TSV files
- Produces distance tables and summary outputs

---

# Recursive Multi-Sample Usage

The user can use a root for a single species that contains multiple individuals:

```
data/
├── 19879801/
│   ├── 801-Cattle.psmc ## Reference
│   ├── 801-Cattle.4.psmc
├── 19999911/
│   ├── 911-Cattle.psmc ## Reference
│   ├── 911-Cattle.4.psmc
```
Therefore, we recommend performing one separate run per each species for an optimal workflow, since each species would need different parameters of **--mu** and **--g**. \nThe reference considered can vary between individuals of the same species, since the each each subfolder must have a reference inside.  

Run:

```bash
python3 PSMC-FAC.py \
  --root data \
  --outdir results_all \
  --base_files 801-Cattle.psmc 911-Cattle.psmc \
  --mu 0.98e-8 \
  --g 5
```

---

# Using a Coverage Table

Instead of `--base_files`, provide:

```
psmc_file                     coverage
sample.psmc                   baseline
sample.4.psmc                 4
sample.5.psmc                 5
```

Run:

```bash
python3 PSMC-FAC.py \
  --root data \
  --outdir results \
  --coverage_table coverage.tsv \
  --mu 0.98e-8 \
  --g 5
```

When `--coverage_table` is used:
- `--base_files` is not required
- Baselines are defined via `baseline` label

---

# Main CLI Arguments

| Argument | Description |
|-----------|------------|
| `--root` | Root directory containing `.psmc` files |
| `--outdir` | Output directory |
| `--base_files` | List of baseline PSMC files |
| `--coverage_table` | TSV mapping files to coverage |
| `--mu` | Mutation rate |
| `--g` | Generation time |
| `--svalue` | PSMC -s parameter (default 100) |
| `--FNR_min` | Minimum FNR to test (default 0) |
| `--FNR_max` | Maximum FNR to test (default 0.99) |
| `--tmin` | Minimum year for custom time vector |
| `--tmax` | Maximum year |
| `--n_timepoints` | Number of timepoints |
| `--distance_scale` | `psmc`, `linear`, or `loglog` |
| `--write_psmc_tsv` | Write corrected TSV per PSMC |
| `--no_curves` | Disable polynomial fitting |
| `--no_trajectory_plots` | Disable trajectory plotting |

---

# Output Files

```
results/
├── full_distance_table.tsv
├── optimal_FNR_summary.tsv
├── polynomial_fits.txt
├── hausdorff_optimal_points.csv
├── frechet_optimal_points.csv
├── trajectories/
└── corrected_tsv/
```

---

# Practical Considerations

- Reference and downsampled files must be in same directory.
- Coverage table must match exactly detected `.psmc` files.
- Polynomial fitting requires ≥ 3 coverage levels.
- Distance scale affects optimal FNR selection.
- Default `psmc` scaling mimics Li & Durbin (2011) plotting.

---

# Interpretation Notes

The optimal FNR minimizes geometric distance but may still show:

- Large SSE in very low coverage samples (<5×)
- Imperfect recovery of reference demographic shape

Inspect both:

- Distance metric
- SSE_log10

for robust interpretation.

---

# Requirements

- Python ≥ 3.8
- numpy
- pandas
- scipy
- matplotlib

Install dependencies:

```bash
pip install numpy pandas scipy matplotlib
```

---

# Citation

If you use PSMC-FAC in your work, please cite:

(Manuscript reference here)

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







+Contact information: 
albanieto@riken.j
francisco.iglesias@univie.ac.at
