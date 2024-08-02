#!/bin/bash

#SBATCH --output=DMRs_all_prediction_ENDTOEND.log                                                    # Standard output and error log.
#SBATCH --qos long                                                              # Partition (queue)
#SBATCH --ntasks=1                                                                      # Run on one mode. 
#SBATCH --cpus-per-task=40                                                              # Number of tasks = cpus. 
#SBATCH --time=4-00:00:00                                                              # Time limit days-hrs:min:sec.
#SBATCH --mem=300gb                                                               # Job memory request.


####### MODULES
source activate group_sRNA

####### VARIABLES
WD="/home/nuezsal/omics_integration_G3/Metiloma"


####### PIPELINE
Rscript 01-DMRcaller_all.R $WD/Results/03-DMRcaller_nuevo $WD/Results/02-Bismark/04-Methylation_extractor $WD/Additional_info/Summary_samples/summary_samples.tsv $SLURM_CPUS_PER_TASK

