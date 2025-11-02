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
├── Run_bwa_clones.sh
├── Run_mpileup_clones.sh
├── Run_bwa_pool.sh
├── Run_mpileup_pool.sh
├── Run_poolFreq.sh
├── singletons_clones.py
├── singletons_summary_positions.py
├── singletons_check_allele_frequency.py
├── Plot_PCA_clones.R
├── Plot_frequencies_from_poolSeq.R
└── README.md
```

---

## Part 1 — Clone Mapping and SNP Discovery

### 1. Mapping individual clones

```bash
sbatch --array=1-30 Run_bwa_clones.sh
```

**Script:** `Run_bwa_clones.sh`  
- Maps each clone against the *D. magna* reference genome (`DmagnaLRV01.fasta`).
- Filters for properly paired, high-quality reads (`-f 0x02 -q 20`).
- Removes PCR duplicates using Picard.
- Output: `*.filtered.sorted.nd.bam` files per sample.

---

### 2. Variant calling

```bash
sbatch Run_mpileup_clones.sh
```

**Script:** `Run_mpileup_clones.sh`  
- Calls SNPs for all 30 clones simultaneously.
- Output: `Daphnia_DmagnaLRV01.vcf.gz`

---

### 3. Filtering SNPs

```bash
vcftools --gzvcf Daphnia_DmagnaLRV01.vcf.gz --recode --remove-indels --min-alleles 2 --minQ 100 --max-missing 1 --stdout > Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf
```

Output: `Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf.gz`

---

### 4. Identify singleton variants

```bash
python singletons_clones.py Daphnia_DmagnaLRV01.mQ100.SNP.NoMissing.vcf > clone_homozygous_singletons_mQ100.txt
```

**Script:** `singletons_clones.py`  
- Parses the VCF file.
- Detects SNPs unique to a single clone (singleton sites).
- Output columns:  
  `Sample  Chrom  Position  Genotype  MajorAllele  UniqueAllele`

---

### 5. Summarize singleton counts

```bash
python singletons_summary_positions.py clone_homozygous_singletons_mQ100.txt
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

### 6. PCA of genotypes

Run PCA in R to visualize population structure.

**Script:** `Plot_PCA_clones.R`  
- Uses the filtered VCF (`*.SNP.NoMissing.vcf.gz`) to perform PCA on genotypes.  
- Distinguishes clones clearly and can reveal hybrid individuals (e.g., DP24 shows intermediate position between DP and DG).

---

## Part 2 — Pool-Seq Analysis

### 1. Map pooled sample(s)

```bash
sbatch Run_bwa_Helene_pool.sh
```

**Script:** `Run_bwa_Helene_pool.sh`  
- Maps pooled sequencing reads to the *D. magna* reference.
- Produces high-quality BAMs with duplicates removed.

---

### 2. Create pileup/VCF with allele depths

```bash
sbatch --array=1-7 Run_mpileup_Helene_pool2.sh
```

**Script:** `Run_mpileup_Helene_pool2.sh`  
- Generates VCFs with per-site allele depth (`AD`) and depth (`DP`) annotations.
- Output: `Daphnia_Helene_DmagnaLRV01_POOL.mpileupRAW_<sample>.vcf`

---

### 3. Estimate clone allele frequencies

```bash
sbatch Run_Helene_poolFreq2.sh
```

**Script:** `Run_Helene_poolFreq2.sh`  
- Iterates over pool VCFs.
- Calls `singletons_check_allele_frequency2.py` to extract allele frequencies of clone-specific singletons.
- Outputs per-sample frequency tables:  
  `Helene_clone_magnapulex_singletons_mQ100_<sample>_ALL_freqs.txt`

---

### 4. Python script details

#### `singletons_check_allele_frequency2.py`
- **Inputs:**
  1. Clone-specific singleton file (`Helene_clone_magnapulex_singletons_mQ100.txt`)
  2. Pool VCF file with AD/DP annotations.
- **Outputs:** Frequency estimates per position and clone.

**Output columns:**
```
Sample  Chromosome  Position  VCF_Genotype  Positions_Genotype
Ref_Allele_Frequency  Alt_Allele_Frequency  Homozygous_or_Heterozygous
```

This allows checking the allele frequency of clone-unique SNPs in the pooled dataset.

---

### 5. Interpretation

- **Sum of mean clone frequencies** ≈ 1.5 → indicates some background noise.  
- **Sum of median frequencies** ≈ 1.0 → confirms accuracy of clone frequency estimates.  
- Deviations likely due to unequal DNA input or differences in individual sizes.

---

## Suggested File Renames (for clarity)

| Current | Suggested |
|----------|------------|
| `Run_bwa_Helene.sh` | `01_map_clones.sh` |
| `Run_mpileup_Helene.sh` | `02_call_variants_clones.sh` |
| `singletons_Helene.py` | `03_extract_singletons.py` |
| `singletons_Helene_summary_positions.py` | `04_summarize_singletons.py` |
| `Run_bwa_Helene_pool.sh` | `05_map_poolseq.sh` |
| `Run_mpileup_Helene_pool2.sh` | `06_call_variants_poolseq.sh` |
| `singletons_check_allele_frequency2.py` | `07_check_allele_frequencies.py` |
| `Run_Helene_poolFreq2.sh` | `08_run_pool_frequency_estimation.sh` |
| `Plot_PCA_Helene.R` | `09_plot_PCA.R` |

---

## Output Summary

| Step | Script | Main Output |
|------|---------|-------------|
| Clone Mapping | `Run_bwa_Helene.sh` | `*.filtered.sorted.nd.bam` |
| Variant Calling | `Run_mpileup_Helene.sh` | `Daphnia_Helene_DmagnaLRV01.vcf.gz` |
| SNP Filtering | `vcftools` | `*.SNP.NoMissing.vcf.gz` |
| Singleton Extraction | `singletons_Helene.py` | `Helene_clone_homozygous_singletons_mQ100.txt` |
| Singleton Summary | `singletons_Helene_summary_positions.py` | Tab summary per clone |
| PCA | `Plot_PCA_Helene.R` | PCA plot |
| Pool Mapping | `Run_bwa_Helene_pool.sh` | `*.filtered.sorted.nd.bam` |
| Pool VCF | `Run_mpileup_Helene_pool2.sh` | `POOL.mpileupRAW_<sample>.vcf` |
| Frequency Estimation | `singletons_check_allele_frequency2.py` | `_ALL_freqs.txt` per sample |

---

## Citation

If you use this pipeline or parts of it, please cite:

> [Your Name] et al., *Daphnia Pool-Seq Pipeline*, GitHub Repository (2025).  
> DOI: [add DOI if archived via Zenodo]

---

## Author

Developed by **Helene [Surname]**  
with pipeline organization and documentation by collaborators.

---

## License

This project is released under the MIT License — feel free to use, modify, and share with attribution.

---

## Notes

For reproducibility, you can archive VCF and summary outputs with metadata describing:
- Reference genome version (`DmagnaLRV01.fasta`)
- Sequencing depth and sample table
- SLURM job IDs for reproducibility

For questions or improvements, open an issue or pull request on GitHub.

---
