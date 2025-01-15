#!/bin/bash

#SBATCH --job-name=more           # Job name to show with squeue
#SBATCH --output=more_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-03:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=50G            # Required Memory per core
#SBATCH --cpus-per-task=4           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   MORE_submit.sh
#
#   This program generate normalize count tables of microRNA, genes and methylome
#   and execute a integration with MLR model of MORE
#
#   Author: Marta Nuñez Salvador
#   Date: 18/12/2024
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate group_sRNA

##################### microRNA ###########################  
# Paths
path_table_sRNA=/home/nuezsal/Omics_integration/sRNA/Results/03-Fusion_count_tables_RF
path_metadata_sRNA=/home/nuezsal/Omics_integration/sRNA/Additional_info
path_out_sRNA=/home/nuezsal/Omics_integration/Integration_MORE/Results/01-sRNA_normalize
path_annot_sRNA=/home/nuezsal/Omics_integration/sRNA/Results/06-miRNAs_grouped_by_family/Group_miRNAs_sig/cume/Omics_project/01-DEA_results_annot

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript normalize_sRNA.R \
            --input $path_table_sRNA \
            --metadata $path_metadata_sRNA \
            --output $path_out_sRNA \
            --annotation $path_annot_sRNA

echo -e "microRNA normalization finish..."





##################### TRANSCRIPTS ###########################            
# Paths
path_table_trans=/home/nuezsal/Omics_integration/RNA_seq/Results/02-Salmon
path_metadata_trans=/home/nuezsal/Omics_integration/RNA_seq/Additional_info
path_out_trans=/home/nuezsal/Omics_integration/Integration_MORE/Results/02-transcripts_normalize
path_annot_trans=/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript normalize_transcripts.R \
            --input $path_table_trans \
            --metadata $path_metadata_trans \
            --output $path_out_trans \
            --annotation $path_annot_trans
            
echo -e "Gene normalization finish..."




##################### METHYLOME ###########################            
# Paths
path_table_met=/home/nuezsal/Omics_integration/Methylome/Results/02-Bismark/04-Methylation_extractor
path_metadata_met=/home/nuezsal/Omics_integration/Methylome/Additional_info
path_out_met=/home/nuezsal/Omics_integration/Integration_MORE/Results/03-methylome_normalize
path_annot_met=/home/nuezsal/Omics_integration/RNA_seq/Results/03-DEA_TH

# Execution 
srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript normalize_methylome.R \
            --input $path_table_met \
            --metadata $path_metadata_met \
            --output $path_out_met \
            --annotation $path_annot_met
            
echo -e "Gene normalization finish..."
