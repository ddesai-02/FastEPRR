#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem 32G
#SBATCH --cpus-per-task=8
#SBATCH --time=11:59:00

module load bcftools

# This command filters for: 
# 1. Biallelic SNPs only
# 2. Minimum Quality 30
# 3. Max missingness 10% (must be in 90% of samples)
# 4. Minor Allele Frequency (MAF) 0.05 (FastEPRR needs polymorphic sites)

bcftools view -v snps -m 2 -M 2 WTD_.vcf.gz | \
bcftools filter -e 'QUAL < 30 || FS > 60.0 || MQ < 40.0' | \
bcftools view -i 'F_MISSING < 0.1 && MAF > 0.05' -Oz -o WTD_filt.vcf.gz

# Index the new file
bcftools index WTD_filt.vcf.gz
