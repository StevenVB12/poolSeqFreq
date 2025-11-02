#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10 
#SBATCH --time=20:00:00 
#SBATCH -A lp_kmma

# Load the programs we will use
module load BWA
module load SAMtools
module load BCFtools

echo "================="

# Sample IDs 
samples=(\
DG12 DG15 DG16 DG29 DG30 \
DG32 DG3 DG33 DG35 DG6 \
DM14 DM18 DM20 DM22 DM25 \
DM27 DM35 DM4 DM5 DM9 \
DP12 DP1 DP18 DP19 DP24 \
DP26 DP33 DP34 DP4 DP8)

# Some folder and file paths to use later
REF=/lustre1/scratch/350/vsc35085/Helene/DmagnaLRV01.fasta
REFNAME=DmagnaLRV01
BWAout=/lustre1/scratch/350/vsc35085/Helene/BAM

# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[*]}
do
ALL_LIST="$ALL_LIST $FILE".$REFNAME.filtered.sorted.nd.bam""
done
eval command=\$$(echo ALL_LIST)

# run mpileup
cd /lustre1/scratch/350/vsc35085/Helene/BAM

bcftools mpileup -O z --threads 10 -f $REF $(echo $command) | bcftools call -m -Oz -o /lustre1/scratch/350/vsc35085/Helene/Daphnia_$REFNAME.vcf.gz 

