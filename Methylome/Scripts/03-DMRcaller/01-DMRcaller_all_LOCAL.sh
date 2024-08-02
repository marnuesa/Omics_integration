#!/bin/bash

#SBATCH --job-name=dmrallL								# Job name.
#SBATCH --output=DMRs_all_prediction_LOCAL.log					# Standard output and error log.
#SBATCH --partition=long								# Partition (queue)
#SBATCH --ntasks=1									# Run on one mode. 
#SBATCH --cpus-per-task=40								# Number of tasks = cpus. 
#SBATCH --time=4-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem-per-cpu=30gb								# Job memory request.


####### MODULES
module load R/4.1.2

####### VARIABLES
WD="/storage/ncRNA/Projects/Integracion_omicas_melon/Metiloma"

####### DIRECTORY
mkdir -p $WD/Results/02-DMRcaller
mkdir -p $WD/Results/02-DMRcaller/LOCAL

####### PIPELINE
Rscript 01-DMRcaller_all.R $WD/Results/02-DMRcaller/LOCAL $WD/Results/01-Bismark/LOCAL/04-Methylation_extractor $WD/Additional_info/Summary_samples/summary_samples.tsv $SLURM_CPUS_PER_TASK
