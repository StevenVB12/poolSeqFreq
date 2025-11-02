#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10 
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA
module load SAMtools/1.9-foss-2018a

echo "================="

# Sample IDs 
samples=(GC159229 GC159230 GC159231 GC159232 GC159233 GC159234)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/lustre1/scratch/350/vsc35085/Helene/DmagnaLRV01.fasta
REFNAME=DmagnaLRV01
BWAout=/lustre1/scratch/350/vsc35085/Helene/BAM

FILE1=/lustre1/scratch/350/vsc35085/data_Helene/$(echo "${samples[ID]}")_combined_R1.fastq.gz
FILE2=/lustre1/scratch/350/vsc35085/data_Helene/$(echo "${samples[ID]}")_combined_R2.fastq.gz

#bwa index $REF

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




