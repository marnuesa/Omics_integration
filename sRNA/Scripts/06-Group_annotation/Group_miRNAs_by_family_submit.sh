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
#   Author: Marta Núñez Salvador
#   Date: 23/09/2024
#   Version: 2.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

############################## With 1 mismatch #################################
# Paths
path_in_dea=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA
path_in_annot=/home/nuezsal/Omics_integration/sRNA/Results/05-Identification_miRNAs
path_out=/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family

# Other variables
data_type=mature

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript Group_miRNAs_by_family.r \
	--dea $path_in_dea \
	--annotation $path_in_annot \
	--type $data_type \
	--output $path_out
    
############################## With 0 mismatch #################################
# Paths
path_in_dea=/home/nuezsal/Omics_integration/sRNA/Results/04-DEA
path_in_annot=/home/nuezsal/Omics_integration/sRNA/Results/05-Identification_miRNAs_simm
path_out=/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family_simm

# Other variables
data_type=mature

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive  Rscript Group_miRNAs_by_family.r \
	--dea $path_in_dea \
	--annotation $path_in_annot \
	--type $data_type \
	--output $path_out

exit 0
