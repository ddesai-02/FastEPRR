import os
import glob
import subprocess
from collections import defaultdict

# --- Configuration (Keep these paths correct) ---
trimmed_dir = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/trimmed"
merge_output_dir = os.path.join(trimmed_dir, "merged_files")

# Ensure the output directory exists
os.makedirs(merge_output_dir, exist_ok=True)
print(f"Trimmed directory: {trimmed_dir}")
print(f"Merge output directory: {merge_output_dir}")
print("-" * 30)

sample_groups = defaultdict(lambda: {'R1': [], 'R2': []})

# --- Step 1: Group Files by Correct Sample ID ---
all_trimmed_files = glob.glob(os.path.join(trimmed_dir, "*_trimmed_R?.fastq.gz"))

for filepath in all_trimmed_files:
    filename = os.path.basename(filepath)
    
    # 1. Determine the Read Group (R1 or R2)
    if '_R1.fastq.gz' in filename:
        read_group = 'R1'
    elif '_R2.fastq.gz' in filename:
        read_group = 'R2'
    else:
        continue 

    # 2. Extract the Sample ID: Everything BEFORE the Lane ID '_L' (e.g., 'Odo_SC1_S6')
    parts_L = filename.split('_L')
    if len(parts_L) < 2:
        print(f"⚠️ Skipping file with no Lane ID ('_L'): {filename}")
        continue
    
    # Isolate the sample_name + sample_index (e.g., 'Odo_SC1_S6')
    sample_plus_index = parts_L[0]
    
    # 3. Final Sample ID: Remove the Sample Index tag ('_SXX')
    # Use rsplit('_S', 1) to split only on the LAST occurrence of '_S'
    final_merge_id = sample_plus_index.rsplit('_S', 1)[0] 
    
    # This logic correctly extracts the unique sample identifier for merging
    # e.g., 'Odo_SC1_S6' -> 'Odo_SC1'
    # e.g., 'Odo_Key3_S3' -> 'Odo_Key3'

    sample_groups[final_merge_id][read_group].append(filepath)

# --- Step 2: Process Samples (Merge or Move) ---
submitted_jobs = 0

for sample_id, files in sample_groups.items():
    # Check if a sample requires merging (i.e., multiple files for R1 or R2)
    if len(files['R1']) > 1 or len(files['R2']) > 1:
        
        # --- Multi-Lane (Merging Required) ---
        print(f"➡️ Multi-lane sample requiring merge: {sample_id} ({len(files['R1'])} lanes)")
        
        # Sort files to ensure R1 and R2 are merged in the same order
        r1_files_sorted = sorted(files['R1'])
        r2_files_sorted = sorted(files['R2'])
        
        # Create the input string for the 'cat' command
        r1_input_str = " ".join(r1_files_sorted)
        r2_input_str = " ".join(r2_files_sorted)
        
        # Define final output names
        final_r1 = os.path.join(merge_output_dir, f"{sample_id}_merged_R1.fastq.gz")
        final_r2 = os.path.join(merge_output_dir, f"{sample_id}_merged_R2.fastq.gz")
        
        # --- Generate and Submit SBATCH script for merging ---
        merge_job_script = f"{sample_id}_merge_job.sh"
        
        script_content = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --job-name=merge_{sample_id}
#SBATCH --output={merge_output_dir}/merge_{sample_id}.out
#SBATCH --error={merge_output_dir}/merge_{sample_id}.err

echo "Starting merge for {sample_id}..."
# Concatenate R1 files
cat {r1_input_str} > {final_r1}

# Concatenate R2 files
cat {r2_input_str} > {final_r2}

if [ $? -eq 0 ]; then
    echo "✅ Merge complete for {sample_id}."
else
    echo "❌ Merge FAILED for {sample_id}."
fi
"""
        with open(merge_job_script, "w") as f:
            f.write(script_content)
        
        try:
            subprocess.run(["sbatch", merge_job_script], check=True, text=True, capture_output=True)
            submitted_jobs += 1
        except subprocess.CalledProcessError as e:
            print(f"❌ Error submitting MERGE job for {sample_id}. SBATCH failed.")

    else:
        # --- Single-Lane (Move/Rename Required) ---
        print(f"✅ Single-lane sample: {sample_id}. No merge needed.")
        
        # Define the target paths for the final files
        final_r1_path = os.path.join(merge_output_dir, f"{sample_id}_merged_R1.fastq.gz")
        final_r2_path = os.path.join(merge_output_dir, f"{sample_id}_merged_R2.fastq.gz")

        # Get the current full paths
        current_r1_path = files['R1'][0]
        current_r2_path = files['R2'][0]

        # Use os.rename to move and rename the files to the final merged directory
        try:
            os.rename(current_r1_path, final_r1_path)
            os.rename(current_r2_path, final_r2_path)
            print(f"   -> Moved and renamed R1/R2 to {merge_output_dir}")
        except OSError as e:
            print(f"❌ Error moving single-lane file for {sample_id}. {e}")


print("-" * 30)
print(f"Total multi-lane samples submitted for merging: {submitted_jobs}")
print("Single-lane files have been moved and renamed.")
print("Wait for SLURM merge jobs to complete before proceeding to alignment.")
