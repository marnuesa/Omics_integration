#!/bin/bash

#SBATCH --job-name=DMRs_plot          # Job name to show with squeue
#SBATCH --output=DMRs_plot_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-01:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=10           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   DMRs_plot_submit.sh
#
#   This program executes the 04-DMRs_plot.R program to perform extract
#   DMRs profile
#
#   Author: Marta Nuñez Salvador
#   Date: 29/08/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_in='/home/nuezsal/Omics_integration/Methylome/Results/03-DMRcaller'
path_out='/home/nuezsal/Omics_integration/Methylome/Results/04-DMRs_plots'

mkdir -p "$path_out"

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript 04-DMRs_plot.R \
            --input $path_in \
            --output $path_out
