#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem 32G
#SBATCH --cpus-per-task=1
#SBATCH --time=11:59:00 
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

#SBATCH --job-name=MergeVCF_L90

echo "Starting GATK MergeVcfs."
echo "Input files found: 34"
echo "Start Time: $(date)"

# --- IMPORTANT: Ensure GATK is loaded or in your PATH ---
# module load gatk/4.x # Example - UNCOMMENT and ADJUST if needed

# 1. MERGEVCFS
gatk --java-options "-Xmx28g -Djava.io.tmpdir=/home/devan/scratch/" MergeVcfs \
    -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069674.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069675.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069676.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069677.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069678.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069679.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069680.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069681.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069682.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069683.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069684.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069685.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069686.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069687.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069688.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069689.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069690.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069691.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069692.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069693.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069694.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069695.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069696.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069697.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069698.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069699.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069700.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069701.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069702.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069703.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069704.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069705.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069706.1.raw.vcf.gz -I /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/final_vcfs/NC_069707.1.raw.vcf.gz \
    -O /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/WTD_.vcf.gz

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Successfully created merged VCF: /home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/WTD_.vcf.gz"
else
    echo "ERROR: GATK MergeVcfs failed. Check stdout/stderr logs."
    exit 1
fi

echo "End Time: $(date)"
