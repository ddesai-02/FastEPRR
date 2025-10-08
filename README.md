# FastEPRR
FastEPRR Workflow from raw fastq files to final recombination rate estimate

STEPS

01_trim.py - trims adaptors, basic QC

02_merge.py - merge individuals run on multiple lanes

03_map.py - maps to reference genome

04_variant_calling.py - calls variants per scaffold per individual

05_gather.py - gathers per individual

06_filter.py - index vcf and filter sites
