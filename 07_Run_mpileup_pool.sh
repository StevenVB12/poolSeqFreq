#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BCFtools/1.15.1-GCC-11.3.0


echo "================="

# Sample IDs 
samples=(GC153502 GC159229 GC159230 GC159231 GC159232 GC159233 GC159234)

# Some folder and file paths to use later
REF=/lustre1/scratch/350/vsc35085/Helene/DmagnaLRV01.fasta
REFNAME=DmagnaLRV01
BWAout=/lustre1/scratch/350/vsc35085/Helene/BAM

# run mpileup
cd /lustre1/scratch/350/vsc35085/Helene/BAM

# need to make sure it outputs AD!!
bcftools mpileup -O v --threads 10 -f $REF $(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam --annotate FORMAT/AD,FORMAT/DP > /lustre1/scratch/350/vsc35085/Helene/Daphnia_Helene_DmagnaLRV01_POOL.mpileupRAW_$(echo "${samples[ID]}").vcf
