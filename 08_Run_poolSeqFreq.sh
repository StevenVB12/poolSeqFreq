#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name freq 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=71:00:00 
#SBATCH --mem-per-cpu=16G
#SBATCH -A lp_svbelleghem

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
#ID=$((SLURM_ARRAY_TASK_ID -1))


#GC153502 GC159229 GC159230 GC159231 GC159232 GC159233 GC159234

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC153502.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC153502_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159229.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159229_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159230.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159230_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159231.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159231_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159232.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159232_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159233.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159233_ALL_freqs.txt

python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
/lustre1/scratch/350/vsc35085/Helene/clone_magnapulicgaleata_singletons_mQ100.txt \
~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159234.vcf > ~/SCRATCH/Helene/clone_magnapulicgaleata_singletons_mQ100_GC159234_ALL_freqs.txt


#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC153502.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC153502_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159229.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159229_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159230.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159230_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159231.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159231_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159232.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159232_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159233.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159233_ALL_freqs.txt

#python /lustre1/scratch/350/vsc35085/data_Helene/singletons_check_allele_frequency.py \
#/lustre1/scratch/350/vsc35085/Helene/clone_pulic_singletons_mQ100.txt \
#~/SCRATCH/Helene/Daphnia_DmagnaLRV01_POOL.mpileupRAW_GC159234.vcf > ~/SCRATCH/Helene/clone_pulic_singletons_mQ100_GC159234_ALL_freqs.txt





