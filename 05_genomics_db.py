import os
import subprocess
import pandas as pd

# --- Configuration ---
# Specified paths
REF_FASTA = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_Ovbor_1.2.fna"
SCAFFOLD_LIST_FILE = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_L90_scaffs.list"
FASTEPRR = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr"
GVCF_DIR = os.path.join(FASTEPRR, "gvcf_files_per_scaffold") # Directory containing input GVCFs

# Output directories
DB_ROOT_DIR = os.path.join(FASTEPRR, "genomics_db_per_scaffold") # Directory for GenomicsDB workspaces
SLURM_SCRIPTS_DIR = os.path.join(FASTEPRR, "slurm_scripts_dbimport")
scratch_dir = "/home/devan/scratch/"

# --- Read Scaffolds into a List ---
try:
    with open(SCAFFOLD_LIST_FILE, 'r') as f:
        # Strip whitespace and ignore empty lines
        SCAFFOLDS = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"ERROR: Scaffold list file not found at {SCAFFOLD_LIST_FILE}")
    exit(1)

# --- Slurm Header ---
MEM_PER_JOB = "110G" 
SLURM_HEADER = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem {MEM_PER_JOB}
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00 
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
"""

# --- Ensure directories exist ---
os.makedirs(DB_ROOT_DIR, exist_ok=True)
os.makedirs(SLURM_SCRIPTS_DIR, exist_ok=True)
print(f"Found {len(SCAFFOLDS)} L90 scaffolds for processing.")
print(f"GenomicsDB workspaces will be created in: {DB_ROOT_DIR}")
print(f"DBImport Slurm scripts will be saved to: {SLURM_SCRIPTS_DIR}")

# --- The Function ---
def combine_scaffold_gvcfs(scaffold, db_root_dir, gvcf_dir, scratch_dir, slurm_scripts_dir):
    """
    Finds all GVCFs for a specific scaffold and creates a GenomicsDBImport Slurm script
    """
    
    # 1. Find all GVCF files for the current scaffold
    gvcf_files = [
        os.path.join(gvcf_dir, f) 
        for f in os.listdir(gvcf_dir) 
        if f.endswith(f"_{scaffold}.g.vcf.gz")
    ]
    
    if not gvcf_files:
        print(f"WARNING: No GVCF files found for scaffold {scaffold}. Skipping.")
        return

    # 2. Build the --variant part of the command
    variant_cmd = ""
    for full_path in gvcf_files:
        variant_cmd += f"--variant {full_path} "
        
    # 3. Define the GenomicsDB output workspace name
    DB_WORKSPACE = os.path.join(db_root_dir, f"db_{scaffold}")
    
    # 4. Construct the  GATK command
    combine_cmd = (
        f"gatk --java-options \"-XX:ParallelGCThreads=1 -Xmx100g -Djava.io.tmpdir={scratch_dir}\" GenomicsDBImport "
        f"{variant_cmd} "
        f"--tmp-dir {scratch_dir} "
        f"--genomicsdb-workspace-path {DB_WORKSPACE} "
        f"-L {scaffold}"
    )
    
    # 5. Define the script name and path
    sh_script_name = f"slurm_dbimport_{scaffold}.sh"
    SH_FILE_PATH = os.path.join(slurm_scripts_dir, sh_script_name)
    
    # --- Write the content to the .sh file ---
    with open(SH_FILE_PATH, 'w') as file:
        file.write(SLURM_HEADER)
        file.write(f'#SBATCH --job-name=DB_{scaffold}\n\n')
        file.write(f'echo "Starting GATK GenomicsDBImport for scaffold: {scaffold}"\n')
        file.write(f'echo "Start Time: $(date)"\n\n')
        
        # Write the main command
        file.write(combine_cmd)
        file.write('\n\n')
        
        # Write the cp command
        file.write(f'cp -a {DB_WORKSPACE} {DB_ROOT_DIR}\n')
        
        file.write('echo "End Time: $(date)"\n')
        file.write(f'echo "Job complete for {scaffold}"\n')
        
    # Make the shell script executable
    os.chmod(SH_FILE_PATH, 0o755)
    print(f"-> Created and made executable: {sh_script_name}")


# --- Run the function for all L90 scaffolds ---
print("\nGenerating GATK GenomicsDBImport scripts per-scaffold...")
for scaffold in SCAFFOLDS:
    combine_scaffold_gvcfs(
        scaffold=scaffold, 
        db_root_dir=DB_ROOT_DIR, 
        gvcf_dir=GVCF_DIR, 
        scratch_dir=scratch_dir, 
        slurm_scripts_dir=SLURM_SCRIPTS_DIR
    )

print("\nScript generation finished for GenomicsDBImport (Step 2).")
print(f"Next step: Submit these scripts from the directory: {SLURM_SCRIPTS_DIR}")
