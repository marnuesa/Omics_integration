#!/bin/bash

#SBATCH --job-name=group          # Job name to show with squeue
#SBATCH --output=group_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-0:10:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=2           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   Group_miRNAs_by_family_submit.sh
#
#   This program executes the Group_miRNAs_by_family.r program.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 12/20/2023
#   Version: 2.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

# Paths
path_in_dea=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA
path_in_annot=/home/nuezsal/Omics_integration/sRNA/Results/05-Identification_miRNAs
path_out=/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family

# Other variables
data_type=mature

# Execution 
Rscript Group_miRNAs_by_family.r \
	--dea $path_in_dea \
	--annotation $path_in_annot \
	--type $data_type \
	--output $path_out
    
exit 0
