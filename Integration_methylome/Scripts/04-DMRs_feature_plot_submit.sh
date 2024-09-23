#!/bin/bash

#SBATCH --job-name=DMRs_features          # Job name to show with squeue
#SBATCH --output=DMRs_features_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-01:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=10           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   DMRs_feature_plot_submit.sh
#
#   This program executes the 04-DMRs_feature_plot.R program to perform extract
#   DMRs features distribution
#
#   Author: Marta Nuñez Salvador
#   Date: 27/08/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_in='/home/nuezsal/Omics_integration/Integration_methylome/Results/02-Total_features_DMRs'
path_out='/home/nuezsal/Omics_integration/Integration_methylome/Results/03-Plots'

mkdir -p "$path_out"

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript 04-DMRs_feature_plot.R \
            --input $path_in \
            --output $path_out
       
