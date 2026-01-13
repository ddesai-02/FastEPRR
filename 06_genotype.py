import os

# --- Configuration ---
REF_FASTA = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_Ovbor_1.2.fna"
SCAFFOLD_LIST_FILE = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr/ref_fasta/WTD_L90_scaffs.list"
FASTEPRR = "/home/devan/projects/def-shaferab/devan/Odocoileus_virginianus/fasteprr"
DB_ROOT_DIR = os.path.join(FASTEPRR, "genomics_db")

# output directories
FINAL_VCF_DIR = os.path.join(FASTEPRR, "final_vcfs") # Directory for final VCFs
SLURM_SCRIPTS_DIR = os.path.join(FASTEPRR, "genotype.log") # Directory for GenotypeGVCFs scripts
scratch_dir = "/home/devan/scratch/" 

# --- Read Scaffolds into a List ---
try:
    with open(SCAFFOLD_LIST_FILE, 'r') as f:
        SCAFFOLDS = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"ERROR: Scaffold list file not found at {SCAFFOLD_LIST_FILE}")
    exit(1)

# --- Slurm Header ---
MEM_PER_JOB = "64G" 
NUM_THREADS = 8
SLURM_HEADER = f"""#!/bin/bash
#SBATCH --account=rrg-shaferab
#SBATCH --mem {MEM_PER_JOB}
#SBATCH --cpus-per-task={NUM_THREADS}
#SBATCH --time=2-23:59:00 
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
"""

# --- Ensure directories exist ---
os.makedirs(FINAL_VCF_DIR, exist_ok=True)
os.makedirs(SLURM_SCRIPTS_DIR, exist_ok=True)
print(f"Genotyping {len(SCAFFOLDS)} scaffolds.")
print(f"Final VCFs will be saved to: {FINAL_VCF_DIR}")
print(f"Genotyping Slurm scripts will be saved to: {SLURM_SCRIPTS_DIR}")


# --- Generate Scripts ---
print("\nGenerating GATK GenotypeGVCFs scripts per-scaffold...")
script_count = 0

# Loop through each scaffold
for scaffold in SCAFFOLDS:
    
    # Define paths
    DB_WORKSPACE = os.path.join(DB_ROOT_DIR, f"db_{scaffold}")
    OUTPUT_VCF = os.path.join(FINAL_VCF_DIR, f"{scaffold}.raw.vcf.gz")
    
    # Script details
    sh_script_name = f"slurm_genotype_{scaffold}.sh"
    SH_FILE_PATH = os.path.join(SLURM_SCRIPTS_DIR, sh_script_name)
    
    # Set Java memory
    java_mem = int(int(MEM_PER_JOB.replace('G', '')) * 0.8)

    # --- Construct the content of the shell script ---
    script_content = f"""{SLURM_HEADER}
#SBATCH --job-name=GTYPE_{scaffold}

echo "Starting GATK GenotypeGVCFs for scaffold: {scaffold}"
echo "Start Time: $(date)"

# --- IMPORTANT: Ensure GATK is loaded or in your PATH ---
# module load gatk/4.x # Example - UNCOMMENT and ADJUST if needed

# 1. GENOTYPEGVCFS (Joint Genotype from the GenomicsDB workspace)
gatk --java-options "-Xmx{java_mem}g -XX:ParallelGCThreads={NUM_THREADS} -Djava.io.tmpdir={scratch_dir}" GenotypeGVCFs \\
    -R {REF_FASTA} \\
    -V gendb://{DB_WORKSPACE} \\
    -O {OUTPUT_VCF} \\
    -L {scaffold}

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Successfully created VCF: {OUTPUT_VCF}"
else
    echo "ERROR: GATK GenotypeGVCFs failed for {scaffold}. Check stdout/stderr logs."
    exit 1
fi

echo "End Time: $(date)"
echo "Job complete for {scaffold}"
"""
    
    # --- Write and submit the .sh file ---
    with open(SH_FILE_PATH, 'w') as f:
        f.write(script_content)
        
    os.chmod(SH_FILE_PATH, 0o755)
    script_count += 1
            
print(f"\nScript generation finished. Total scripts created: {script_count} (one per scaffold).")
print(f"Next step: Submit these scripts from the directory: {SLURM_SCRIPTS_DIR}")
