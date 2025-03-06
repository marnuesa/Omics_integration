#!/bin/bash

#SBATCH --job-name=Integration_mirnas.log           # Job name to show with squeue
#SBATCH --output=Integration_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-03:00:00          # Time limit to execute the job
#SBATCH --mem=50G            # Required Memory per core
#SBATCH --cpus-per-task=24           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   Integration_miRNAs.sh
#
#   This program executes the miRNA-targets_profile.r program to study the 
#   expression profile of each microRNA and their targets
#
#   Author: Marta Nuñez Salvador
#   Date: 26/07/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_micro="/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family"
path_transcripts="/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA/01-DEA_raw/cume/Omics_project"
path_out="/home/nuezsal/Omics_integration/Integration_microRNA/Results"
path_ai="/home/nuezsal/Omics_integration/Integration_microRNA/Additional_info"


# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript miRNA-targets_profile.R \
            --microrna $path_micro \
            --transcripts $path_transcripts \
            --output $path_out \
            --additional $path_ai 
