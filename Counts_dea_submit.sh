#!/bin/bash

#SBATCH --job-name=Counts_dea   # Job name to show with squeue
#SBATCH --output=Counts_dea_%j.out  # Output file
#SBATCH --ntasks=1                   # Maximum number of cores to use
#SBATCH --time=00-07:00:00           # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=2            # CPUs assigned per task.
#SBATCH --qos=short                  # QoS: short, medium, long, long-mem

#******************************************************************************
#
# Counts_dea_submit.sh
#
# This script executes the Counts_dea.R script to process omics data
# (methylome, microRNA, and RNA-seq) and generate differential expression counts.
#
# Author: Marta Núñez Salvador
# Date: 22/09/2024
# Version: 1.0
#
#******************************************************************************

# Load necessary modules
module load anaconda
conda activate group_sRNA

# Paths
methylome_path="/home/nuezsal/Omics_integration/Methylome/Results/04-DMRs_plots"
microRNA_path="/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family/Group_miRNAs_sig/cume/Omics_project/01-DEA_results_annot"
rnaseq_path="/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH/02-DEA_sig/cume/Omics_project"
output_path="/home/nuezsal/Omics_integration/Summary"

# Execution
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript Counts_dea.R \
  --methylome_path $methylome_path \
  --microrna_path $microRNA_path \
  --rnaseq_path $rnaseq_path \
  --output_path $output_path
