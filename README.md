# Pool-Seq-Freq Pipeline

This repository contains scripts and instructions for processing Pool-Seq and individual clone sequencing data for *Daphnia magna* and related species.  
The workflow includes read mapping, variant calling, singleton detection, and allele frequency estimation in pooled samples. While the pipeline here is specific to *Daphnia magna*, it could be easily modified for other clonal species.

---

## Overview

The pipeline is divided into two main parts:

1. **Part 1 — Individual Clone Processing**
   - Map 30 individual clone sequencing reads against the *D. magna* reference genome.
   - Call SNPs using `bcftools mpileup` and filter high-quality, non-missing positions.
   - Identify clone-specific singleton SNPs (unique to one clone).
   - Summarize singleton counts and visualize population structure via PCA.

2. **Part 2 — Pool-Seq Processing**
   - Map pooled reads against the *D. magna* reference genome.
   - Generate allele frequency information using `bcftools mpileup` (annotating AD and DP fields).
   - Estimate clone frequencies in the pool using singleton markers.

---

## Dependencies

The scripts assume an HPC environment using SLURM and the following modules/tools:

| Tool | Version | Purpose |
|------|----------|----------|
| BWA | ≥0.7.17 | Mapping reads |
| SAMtools | ≥1.9 | BAM filtering/sorting |
| Picard | ≥2.18 | Duplicate removal |
| BCFtools | ≥1.15 | Variant calling |
| Python 3 | ≥3.8 | Data parsing and filtering |
| R | ≥4.0 | PCA and visualization |

---

## Directory Structure

```
.
├── 01_Run_bwa_clones.sh
├── 02_Run_mpileup_clones.sh
├── 03_singletons_clones.py
├── 04_singletons_summary_positions.py
├── 05_Plot_PCA_clones.R
├── 06_Run_bwa_pool.sh
├── 07_Run_mpileup_pool.sh
├── 08_Run_poolSeqFreq.sh
├── singletons_check_allele_frequency.py
├── 09_Plot_frequencies_from_poolSeq.R
└── README.md
```

---

## Part 1 — Clone Mapping and SNP Discovery

### 1. Mapping individual clones

```bash
sbatch --array=1-30 01_Run_bwa_clones.sh
```

**Script:** `01_Run_bwa_clones.sh`  
- Maps each clone against the *D. magna* reference genome (`DmagnaLRV01.fasta`).
- Filters for properly paired, high-quality reads (`-f 0x02 -q 20`).
- Removes PCR duplicates using Picard.
- Output: `*.filtered.sorted.nd.bam` files per sample.

---

### 2.1 Variant calling

```bash
sbatch 02_Run_mpileup_clones.sh
```

**Script:** `02_Run_mpileup_clones.sh`  
- Calls SNPs for all 30 clones simultaneously.
- Output: `Daphnia_DmagnaLRV01.vcf.gz`

---

### 2.2 Filtering SNPs

```bash
vcftools --gzvcf Daphnia_DmagnaLRV01.vcf.gz --recode --remove-indels --min-alleles 2 --minQ 100 --max-missing 1 --stdout > Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf
```

Output: `Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf.gz`

---

### 3. Identify singleton variants

```bash
python 03_singletons_clones.py Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf > clone_homozygous_singletons_mQ100.txt
```

**Script:** `03_singletons_clones.py`  
- Parses the VCF file.
- Detects SNPs unique to a single clone (singleton sites).
- Output columns:  
  `Sample  Chrom  Position  Genotype  MajorAllele  UniqueAllele`

---

### 4. Summarize singleton counts

```bash
python 04_singletons_summary_positions.py clone_homozygous_singletons_mQ100.txt
```

Output example:

```
DG30    178
DG29    223
DM4     303
...
```

Each number represents the count of unique singleton positions per clone.

---

### 5. PCA of genotypes

Run PCA in R to visualize population structure.

**Script:** `05_Plot_PCA_clones.R`  
- Uses the filtered VCF (`*.SNP.NoMissing.vcf.gz`) to perform PCA on genotypes.  
- Distinguishes clones clearly and can reveal hybrid individuals (e.g., DP24 shows intermediate position between DP and DG).

---



## Part 2 — Pool-Seq Analysis

### 1. Map pooled sample(s)

```bash
sbatch 06_Run_bwa_pool.sh
```

**Script:** `06_Run_bwa_pool.sh`  
- Maps pooled sequencing reads to the *D. magna* reference.
- Produces high-quality BAMs with duplicates removed.

---

### 2. Create pileup/VCF with allele depths

```bash
sbatch --array=1-7 07_Run_mpileup_pool.sh
```

**Script:** `07_Run_mpileup_pool.sh`  
- Generates VCFs with per-site allele depth (`AD`) and depth (`DP`) annotations.
- Output: `Daphnia_DmagnaLRV01_POOL.mpileupRAW_<sample>.vcf`

---

### 3. Estimate clone allele frequencies

```bash
sbatch 08_Run_poolSeqFreq.sh
```

**Script:** `08_Run_poolSeqFreq.sh`  
- Iterates over pool VCFs.
- Calls `singletons_check_allele_frequency.py` to extract allele frequencies of clone-specific singletons.
- Outputs per-sample frequency tables:  
  `clone_magnapulicgaleata_singletons_mQ100_<sample>_ALL_freqs.txt`

---

### 4. Python script details

#### `singletons_check_allele_frequency.py`
- **Inputs:**
  1. Clone-specific singleton file (`clone_magnapulicgaleata_singletons_mQ100.txt`)
  2. Pool VCF file with AD/DP annotations.
- **Outputs:** Frequency estimates per position and clone.

**Output columns:**
```
Sample  Chromosome  Position  VCF_Genotype  Positions_Genotype
Ref_Allele_Frequency  Alt_Allele_Frequency  Homozygous_or_Heterozygous
```

This allows checking the allele frequency of clone-unique SNPs in the pooled dataset.

---

### 5. Calculate and Visualize Clone Frequencies

```bash
Rscript 09_Calculate_frequencies_Helene.R
```

**Script:** `09_Calculate_frequencies_Helene.R`  
- Aggregates all per-sample frequency tables (`*_ALL_freqs.txt`) into a single dataset.  
- Corrects allele frequencies for zygosity (`Alt_Allele_Frequency_corrected`).  
- Produces histograms of corrected allele frequencies for each pooled sample.  
- Calculates mean, median, and 95% confidence intervals per sample.  
- **Outputs:**  
  - `summary_stats.txt` — table of clone-level summary statistics  
  - Frequency distribution plots per sample  

---

### 6. Interpretation

- **Sum of median frequencies** ≈ 1.0 → confirms accuracy of clone frequency estimates.  
- Deviations likely due to unequal DNA input or differences in individual sizes.


---

## Citation

If you use this pipeline or parts of it, please cite:

Vanvelk, H., Govaert, L., Van Belleghem, S., Matthews, B., Spaak, P., & De Meester, L. (2025). Genetic diversity alters zooplankton community assembly responses to fish predation.

---

## Author

Developed by **Héléne Vanvelk and Steven Van Belleghem**  

---

## License

This project is released under the MIT License — feel free to use, modify, and share with attribution.

---

## Notes

Daphnia magna reference genome used: https://pmc.ncbi.nlm.nih.gov/articles/PMC10570034/ 

---
