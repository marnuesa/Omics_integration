#!/bin/bash

#SBATCH --job-name=go_enrichment      # Job name to show with squeue
#SBATCH --output=go_enrichment_%j.out # Output file
#SBATCH --ntasks=1                    # Maximum number of cores to use
#SBATCH --time=00-03:00:00             # Time limit to execute the job
#SBATCH --mem-per-cpu=50G              # Required Memory per core
#SBATCH --cpus-per-task=2              # CPUs assigned per task.
#SBATCH --qos=short                    # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   GO_Enrichment_submit.sh
#
#   This script executes the GO_Enrichment_Analysis.R script to perform Gene 
#   Ontology enrichment analysis on differentially expressed genes.
#
#   Author: Your Name
#   Date: 19/02/2025
#   Version: 1.0 
#
#******************************************************************************

# Load necessary modules
module load anaconda
source activate GO_enrich

# Paths
transcripts_path="/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH/02-DEA_sig/cume/Omics_project"
annotation_file="/home/nuezsal/Omics_integration/RNA_seq/Additional_info/DHL92_gene_description_v4.txt"
output_path="/home/nuezsal/Omics_integration/RNA_seq/Results/Enrich_GO"

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript Go_Enrichment_Analysis.R \
            --input $transcripts_path \
            --annotation $annotation_file \
            --output $output_path

