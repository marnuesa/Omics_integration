#!/bin/bash

#SBATCH --job-name=more # Job name to show with squeue
#SBATCH --output=more_%j.out      # Output file
#SBATCH --ntasks=1                              # Maximum number of cores to use
#SBATCH --time=00-07:00:00                      # Time limit to execute the job
#SBATCH --mem-per-cpu=50G                       # Required Memory per core
#SBATCH --cpus-per-task=2                       # CPUs assigned per task.
#SBATCH --qos=short                             # QoS: short,medium,long,long-mem

#******************************************************************************
#
# MORE_submit.sh
#
# This script executes the Methylome_MicroRNA_Associations.R script to analyze
# relationships between methylation, microRNA, and transcriptomic data.
#
# Author: Marta Núñez Salvador
# Date: 19/02/2025
# Version: 1.0
#
#******************************************************************************

# Load necessary modules
module load anaconda
source activate group_sRNA

# Paths
metadata_path="/home/nuezsal/Omics_integration/Integration_MORE/Additional_info"
methylome_associations_file="/home/nuezsal/Omics_integration/Integration_MORE/Results/Bed_files/CMelon_DHL92_v4_upstream_filter_uniq.bed"
microRNA_associations_file="/home/nuezsal/Omics_integration/Integration_microRNA/Additional_info/Targets_orthologues.txt"
input_transcripts_file="/home/nuezsal/Omics_integration/Integration_MORE/Results/02-transcripts_normalize/01-Tables/Genes_normalize_counts.tsv"
input_microRNA_file="/home/nuezsal/Omics_integration/Integration_MORE/Results/01-sRNA_normalize/01-Tables/sRNA_normalize_counts.tsv"
input_methylome_file="/home/nuezsal/Omics_integration/Integration_MORE/Results/03-methylome_normalize/01-Tables/methylation_normalize_counts_upstream.tsv"
DEG_path="/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA/02-DEA_sig/cume/Omics_project"
output_path="/home/nuezsal/Omics_integration/Integration_MORE/Results/04-MORE"

# Execution
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript MORE_integration.R \
  --metadata $metadata_path \
  --methylome_associations $methylome_associations_file \
  --microRNA_associations $microRNA_associations_file \
  --input_transcripts $input_transcripts_file \
  --input_microRNA $input_microRNA_file \
  --input_methylome $input_methylome_file \
  --input_DEG $DEG_path \
  --output $output_path
