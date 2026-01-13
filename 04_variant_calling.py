import os
import glob
import subprocess
import re

# --- Configuration ---

project_root = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/trimmed"
gatk_base_dir = os.path.join(project_root, "gatk_variants")
fixed_bam_dir = os.path.join(gatk_base_dir, "fixed_bams")
gvcf_dir = os.path.join(gatk_base_dir, "gvcf_files")
gvcf_split_dir = os.path.join(gvcf_dir, "split_scaffolds") 
REFERENCE_FASTA = os.path.join(os.path.dirname(project_root), "ref_fasta/WTD_Ovbor_1.2.fna")

# List of scaffolds
SCAFFOLD_LIST = [
    "NC_069674.1", "NC_069675.1", "NC_069676.1", "NC_069677.1", 
    "NC_069678.1", "NC_069679.1", "NC_069680.1", "NC_069681.1", 
    "NC_069682.1", "NC_069683.1", "NC_069684.1", "NC_069685.1", 
    "NC_069686.1", "NC_069687.1", "NC_069688.1", "NC_069689.1", 
    "NC_069690.1", "NC_069691.1", "NC_069692.1", "NC_069693.1", 
    "NC_069694.1", "NC_069695.1", "NC_069696.1", "NC_069697.1", 
    "NC_069698.1", "NC_069699.1", "NC_069700.1", "NC_069701.1", 
    "NC_069702.1", "NC_069703.1", "NC_069704.1", "NC_069705.1", 
    "NC_069706.1", "NC_069707.1", "NC_069708.1", "NC_069709.1"
]

# Ensure output directories exist
os.makedirs(gvcf_split_dir, exist_ok=True)

print("--- Submitting HaplotypeCaller Jobs ---")

# Identify sample IDs from BAM files
fixed_bam_files = glob.glob(os.path.join(fixed_bam_dir, "*_fixed_rg.bam"))
sample_bams = {os.path.basename(f).replace("_fixed_rg.bam", ""): f for f in fixed_bam_files}
sample_ids = list(sample_bams.keys())

if not sample_ids:
    print(f"❌ ERROR: No fixed BAM files found in {fixed_bam_dir}.")
    exit(1)

total_jobs = len(sample_ids) * len(SCAFFOLD_LIST)
print(f"Submitting {total_jobs} total HaplotypeCaller jobs...")

submitted_jobs = []

for sample_id in sample_ids:
    bam_path = sample_bams[sample_id]
    
    for scaffold in SCAFFOLD_LIST:
        gvcf_output = os.path.join(gvcf_split_dir, f"{sample_id}_{scaffold}.g.vcf.gz")
        job_script_path = os.path.join(gatk_base_dir, f"{sample_id}_{scaffold}_02_HC.sh")
        
        # SBATCH script content
        script_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --job-name=HC_{sample_id}_{scaffold}
#SBATCH --output={gvcf_split_dir}/{sample_id}_{scaffold}.out
#SBATCH --error={gvcf_split_dir}/{sample_id}_{scaffold}.err

module load StdEnv/2020
module load java/14.0.2
module load gatk/4.2.5.0

echo "Starting HaplotypeCaller for {sample_id} on {scaffold}..."

gatk HaplotypeCaller \\
    -R {REFERENCE_FASTA} \\
    -I {bam_path} \\
    -O {gvcf_output} \\
    -L {scaffold} \\
    --emit-ref-confidence GVCF \\
    --max-alternate-alleles 3 \\
    -ERC GVCF

if [ $? -ne 0 ]; then
    echo "HaplotypeCaller FAILED for {sample_id} on {scaffold}."
    exit 1
fi

echo "✅ HaplotypeCaller complete for {sample_id} on {scaffold}."
"""

        with open(job_script_path, "w") as f:
            f.write(script_content)
        
        # Submit the job
        try:
            result = subprocess.run(["sbatch", job_script_path], check=True, text=True, capture_output=True)
            job_id_match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if job_id_match:
                submitted_jobs.append(job_id_match.group(1))
        except subprocess.CalledProcessError as e:
            print(f"❌ Error submitting job for {sample_id} on {scaffold}. SBATCH failed. Error: {e.stderr}")

print(f"\nTotal scatter jobs submitted: {len(submitted_jobs)}.")
