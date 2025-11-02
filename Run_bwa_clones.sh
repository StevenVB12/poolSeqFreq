#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10 
#SBATCH --time=8:00:00 
#SBATCH -A lp_kmma

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA
module load SAMtools

echo "================="

# Sample IDs 
samples=(\
DG12 DG15 DG16 DG29 DG30 \
DG32 DG3 DG33 DG35 DG6 \
DM14 DM18 DM20 DM22 DM25 \
DM27 DM35 DM4 DM5 DM9 \
DP12 DP1 DP18 DP19 DP24 \
DP26 DP33 DP34 DP4 DP8)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/lustre1/scratch/350/vsc35085/Helene/DmagnaLRV01.fasta
REFNAME=DmagnaLRV01
BWAout=/lustre1/scratch/350/vsc35085/Helene/BAM

FILE1=/lustre1/scratch/350/vsc35085/Helene/Genomics/$(echo "${samples[ID]}")_1.fq.gz
FILE2=/lustre1/scratch/350/vsc35085/Helene/Genomics/$(echo "${samples[ID]}")_2.fq.gz

# Run BWA mapping
bwa mem -t 10 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove PCR duplicates
java -jar /vsc-hard-mounts/leuven-data/350/vsc35085/programs/picard.jar MarkDuplicates \
-I $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
-O $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-M $BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt \
-ASSUME_SORTED true

# Remove intermediate files
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

echo "================="




