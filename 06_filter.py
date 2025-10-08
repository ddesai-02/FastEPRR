import os
import glob
import subprocess
import re

# ======================================================================
# --- CONFIGURATION ---
# IMPORTANT: Must match the paths in previous scripts
# ======================================================================
project_root = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/trimmed"
gatk_base_dir = os.path.join(project_root, "gatk_variants")
gvcf_dir = os.path.join(gatk_base_dir, "gvcf_files")
gvcf_merged_dir = os.path.join(gvcf_dir, "merged_samples") 
vcf_filtered_dir = os.path.join(gatk_base_dir, "filtered_vcf")
REFERENCE_FASTA = os.path.join(os.path.dirname(project_root), "ref_fasta/WTD_Ovbor_1.2.fna")

# Ensure output directory exists
os.makedirs(vcf_filtered_dir, exist_ok=True)

# ======================================================================
# --- WORKFLOW ---
# ======================================================================
print("--- 4. Submitting GVCF Indexing and VCF Filtering Jobs ---")

# 1. Identify all merged GVCF files (from successful Gather step)
merged_gvcf_files = glob.glob(os.path.join(gvcf_merged_dir, "*.g.vcf.gz"))

if not merged_gvcf_files:
    print(f"❌ ERROR: No merged GVCF files found in {gvcf_merged_dir}. Run 03_submit_gather_gvcf.py first.")
    exit(1)

index_jobs_info = []

# --- Stage 4a: Index Merged GVCFs (MANDATORY FIX) ---
print("Submitting Indexing jobs for merged GVCFs...")
for gvcf_path in merged_gvcf_files:
    sample_id = os.path.basename(gvcf_path).replace(".g.vcf.gz", "")
    index_job_script = os.path.join(gatk_base_dir, f"{sample_id}_04a_index_gvcf.sh")
    
    index_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --job-name=Index_GVCF_{sample_id}
#SBATCH --output={vcf_filtered_dir}/{sample_id}_index.out
#SBATCH --error={vcf_filtered_dir}/{sample_id}_index.err

module load StdEnv/2020
module load bcftools/1.11

echo "Starting BCFtools Index for {sample_id}..."

# Use BCFtools index to create the .tbi file required by GATK SelectVariants
bcftools index -t {gvcf_path}

if [ $? -eq 0 ]; then
    echo "✅ Indexing complete for {sample_id}."
else
    echo "❌ Indexing FAILED for {sample_id}."
    exit 1
fi
"""
    with open(index_job_script, "w") as f:
        f.write(index_content)
    
    # Submit the Index job and capture the Job ID
    try:
        result = subprocess.run(["sbatch", index_job_script], check=True, text=True, capture_output=True)
        job_id_match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if job_id_match:
            job_id = job_id_match.group(1)
            index_jobs_info.append({'id': job_id, 'sample': sample_id, 'gvcf_path': gvcf_path})
            print(f"   -> Submitted Index job for {sample_id} (ID: {job_id})")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error submitting Index job for {sample_id}. SBATCH failed. Error: {e.stderr}")

print("\n--- Submitting Chained Select/Filter Jobs ---")
# --- Stage 4b: Select Variants and Filter (Chained) ---
submitted_filter_jobs = []

for job_info in index_jobs_info:
    sample_id = job_info['sample']
    index_job_id = job_info['id']
    gvcf_path = job_info['gvcf_path']
    vcf_output_path = os.path.join(vcf_filtered_dir, f"{sample_id}_variants.vcf.gz")
    
    vcf_job_script = os.path.join(gatk_base_dir, f"{sample_id}_04b_vcf_filter.sh")
    
    # Filtering Criteria: QUAL>30 (variant quality), FORMAT/DP>5 (min depth), MQ>40 (mapping quality)
    filter_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=VCF_Filter_{sample_id}
#SBATCH --output={vcf_filtered_dir}/{sample_id}.out
#SBATCH --error={vcf_filtered_dir}/{sample_id}.err
#SBATCH --dependency=afterok:{index_job_id}

module load StdEnv/2020
module load java/14.0.2
module load gatk/4.2.5.0
module load bcftools/1.11

echo "Starting SelectVariants and Filtering for {sample_id}..."

# 1. GATK SelectVariants: Extracts only variant sites (requires indexed GVCF)
# 2. BCFtools Filter: Applies sample-level quality filters
# 3. BCFtools Index: Indexes the final VCF
gatk SelectVariants \\
    -R {REFERENCE_FASTA} \\
    -V {gvcf_path} \\
    -O /dev/stdout \\
    --select-type-to-include SNP \\
    --select-type-to-include INDEL \\
    | bcftools filter \\
    -i 'QUAL>30 && FORMAT/DP>5 && MQ>40' \\
    -O z \\
    -o {vcf_output_path}

if [ $? -eq 0 ]; then
    echo "✅ VCF filtering complete for {sample_id}. Indexing final VCF..."
    bcftools index {vcf_output_path}
    if [ $? -eq 0 ]; then
        echo "✅ Final VCF Indexing complete."
    else
        echo "❌ Final VCF Indexing FAILED."
        exit 1
    fi
else
    echo "❌ VCF Conversion/Filtering FAILED for {sample_id}. Check log file."
    exit 1
fi
"""
    with open(vcf_job_script, "w") as f:
        f.write(filter_content)
    
    # Submit the VCF job
    try:
        result = subprocess.run(["sbatch", vcf_job_script], check=True, text=True, capture_output=True)
        job_id_match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if job_id_match:
            submitted_filter_jobs.append(job_id_match.group(1))
            print(f"   -> Submitted Filter job for {sample_id} (ID: {job_id_match.group(1)}), depends on Index job {index_job_id}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error submitting VCF Filter job for {sample_id}. SBATCH failed. Error: {e.stderr}")

print("-" * 30)
print(f"Total Index jobs submitted: {len(index_jobs_info)}")
print(f"Total Filter jobs submitted: {len(submitted_filter_jobs)}")
print("The Select/Filter jobs are chained to run *only* after their corresponding Index job finishes. You now have a complete, segmented workflow!")
