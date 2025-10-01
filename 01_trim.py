import os
import glob
import subprocess

# --- Configuration ---
# Directory where your FASTQ files are located
fastq_dir = "/home/devan/projects/rrg-shaferab/DEER_FASTQ"
# Output directory for trimmed FASTQ files
output_dir = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/trimmed"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)
print(f"FASTQ directory: {fastq_dir}")
print(f"Output directory: {output_dir}")
print("-" * 30)

# Find all R1 files matching the actual naming pattern
r1_files = glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz"))

if not r1_files:
    # This check should now pass if the files are there
    print(f"Error: No *_R1_001.fastq.gz files found in {fastq_dir}. Check the path/pattern.")

for r1 in r1_files:
    # Extract the sample name by removing the correct suffix
    # The sample name will be everything up to "_R1_001.fastq.gz"
    sample = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
    
    # Construct the full path for R2
    r2 = os.path.join(fastq_dir, f"{sample}_R2_001.fastq.gz")

    if not os.path.exists(r2):
        print(f"⚠️ Missing R2 file: {r2}, skipping sample {sample}.")
        continue

    # Define job script filename
    job_script_path = f"{sample}_fastp.sh"
    
    # Define full paths for output files
    out_r1 = os.path.join(output_dir, f"{sample}_trimmed_R1.fastq.gz")
    out_r2 = os.path.join(output_dir, f"{sample}_trimmed_R2.fastq.gz")
    html_report = f"{sample}_fastp.html" 

    # --- Generate the SBATCH script ---
    script_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:59:00
#SBATCH --job-name=fastp_{sample}
#SBATCH --output=fastp_{sample}.out
#SBATCH --error=fastp_{sample}.err

# Load required modules
module load fastp

# Step 1: Quality control and trimming (All paths are absolute from the Python script)
fastp -i {r1} -I {r2} -o {out_r1} -O {out_r2} --detect_adapter_for_pe --thread 1 --html {html_report}
"""

    with open(job_script_path, "w") as f:
        f.write(script_content)

    # --- Submit the job ---
    try:
        # Use shell=True for simple commands in environments like this, 
        # though subprocess.run(["sbatch", job_script_path], ...) is generally safer.
        # Sticking with the list format as it's better practice:
        subprocess.run(["sbatch", job_script_path], check=True, text=True, capture_output=True)
        print(f"✅ Submitted job for: {sample}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error submitting job for {sample}. SBATCH failed.")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")

print("-" * 30)
print("Finished generating and submitting all jobs.")
