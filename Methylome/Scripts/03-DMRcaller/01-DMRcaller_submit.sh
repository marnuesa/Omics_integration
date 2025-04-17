#!/bin/bash
#SBATCH --output=DMRs_%j.log                                                    # Standard output and error log.
#SBATCH --qos long                                                              # Partition (queue)
#SBATCH --ntasks=1                                                                      # Run on one mode. 
#SBATCH --cpus-per-task=40                                                              # Number of tasks = cpus. 
#SBATCH --time=10-00:00:00                                                              # Time limit days-hrs:min:sec.
#SBATCH --mem=300gb                                                               # Job memory request.

################################################################################
#                       03-DMRcaller.sh
#
#   This script runs DMRcaller using an R script to identify differentially
#   methylated regions (DMRs) from bisulfite sequencing data.
#
#   Input  : Methylation extractor output, sample summary table
#   Output : DMR tables and visualizations
#
################################################################################

####### MODULES
module load anaconda
source activate group_sRNA

####### VARIABLES
path="/home/nuezsal/Omics_integration/Methylome"


####### PIPELINE
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript 01-DMRcaller.R \
	-o $path/Results/03-DMRcaller/ \
	-i $path/Results/02-Bismark/04-Methylation_extractor \
	-s $path/Additional_info/Summary_samples/summary_samples.tsv \
	-c $SLURM_CPUS_PER_TASK
