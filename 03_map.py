import os
import glob
import subprocess
import re

# ======================================================================
# --- CONFIGURATION ---
# IMPORTANT: Adjust these paths to match your project structure.
# ======================================================================
project_root = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/trimmed"
# Directory containing R1/R2 FASTQ files (e.g., /path/to/trimmed_reads)
fastq_dir = os.path.join(project_root, "trimmed_reads") 
# Directory where BAM files will be stored
mapping_dir = os.path.join(project_root, "bwa_mapping")
# Reference FASTA file
REFERENCE_FASTA = os.path.join(os.path.dirname(project_root), "ref_fasta/WTD_Ovbor_1.2.fna")
# New directory for GATK outputs (where fixed BAMs will land)
gatk_base_dir = os.path.join(project_root, "gatk_variants")
fixed_bam_dir = os.path.join(gatk_base_dir, "fixed_bams")

# Ensure all output directories exist
os.makedirs(mapping_dir, exist_ok=True)
os.makedirs(fixed_bam_dir, exist_ok=True)

# ======================================================================
# --- WORKFLOW ---
# ======================================================================
print("--- 1. Submitting BWA Mapping and BAM Fixing Jobs ---")

# 1. Locate FASTQ pairs (Assumes files are named SAMPLE_R1.fastq.gz and SAMPLE_R2.fastq.gz)
fastq_r1_files = glob.glob(os.path.join(fastq_dir, "*_R1.fastq.gz"))
sample_ids = [os.path.basename(f).replace("_R1.fastq.gz", "") for f in fastq_r1_files]

if not sample_ids:
    print(f"❌ ERROR: No FASTQ R1 files found in {fastq_dir}. Check the path and naming.")
    exit(1)

print(f"Found {len(sample_ids)} samples to process: {', '.join(sample_ids[:3])}...")

submitted_jobs = []

for sample_id in sample_ids:
    fq1 = os.path.join(fastq_dir, f"{sample_id}_R1.fastq.gz")
    fq2 = os.path.join(fastq_dir, f"{sample_id}_R2.fastq.gz")
    
    # Paths for intermediate and final files
    sam_output = os.path.join(mapping_dir, f"{sample_id}.sam")
    sorted_bam = os.path.join(mapping_dir, f"{sample_id}_sorted.bam")
    fixed_bam = os.path.join(fixed_bam_dir, f"{sample_id}_fixed_rg.bam")
    
    job_script_path = os.path.join(gatk_base_dir, f"{sample_id}_01_map_fix.sh")
    
    # SBATCH script content
    script_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --job-name=MapFix_{sample_id}
#SBATCH --output={mapping_dir}/{sample_id}.out
#SBATCH --error={mapping_dir}/{sample_id}.err

module load StdEnv/2020
module load bwa/0.7.17
module load samtools/1.17
module load java/14.0.2
module load gatk/4.2.5.0 # Use consistent GATK version

echo "Starting BWA Mapping for {sample_id}..."

# 1. BWA Mapping
# Assumes the reference FASTA has been indexed (bwa index)
bwa mem -t 8 -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' \\
    {REFERENCE_FASTA} {fq1} {fq2} > {sam_output}
if [ $? -ne 0 ]; then echo "BWA Mapping FAILED."; exit 1; fi

echo "Sorting and converting to BAM..."
# 2. SAMtools Sort
samtools view -@ 8 -bS {sam_output} \\
    | samtools sort -@ 8 -o {sorted_bam}
if [ $? -ne 0 ]; then echo "SAMtools Sort FAILED."; exit 1; fi

# 3. GATK AddOrReplaceReadGroups (to ensure proper read group headers for HaplotypeCaller)
echo "Adding/Replacing Read Groups and indexing BAM..."
gatk AddOrReplaceReadGroups \\
    -I {sorted_bam} \\
    -O {fixed_bam} \\
    -RGID {sample_id} -RGLB {sample_id}_lib -RGPL ILLUMINA -RGSM {sample_id} -RGPU unit1 \\
    --CREATE-INDEX true
if [ $? -ne 0 ]; then echo "GATK ReadGroups FAILED."; exit 1; fi

# 4. Clean up intermediate files
rm {sam_output}
rm {sorted_bam}

echo "✅ Mapping and BAM fixing complete for {sample_id}."
"""
    
    with open(job_script_path, "w") as f:
        f.write(script_content)
    
    # Submit the job
    try:
        result = subprocess.run(["sbatch", job_script_path], check=True, text=True, capture_output=True)
        job_id_match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if job_id_match:
            submitted_jobs.append(job_id_match.group(1))
            print(f"   -> Submitted {sample_id} (Job ID: {job_id_match.group(1)})")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error submitting job for {sample_id}. SBATCH failed. Error: {e.stderr}")

print(f"\nTotal jobs submitted: {len(submitted_jobs)}. Check your queue.")
