#!/bin/bash
#SBATCH --output=bismark_%j.log							# Standard output and error log.
#SBATCH --qos=medium								# Partition (queue)
#SBATCH --ntasks=4									# Run on one mode. 
#SBATCH --cpus-per-task=26								# Number of tasks = cpus. 
#SBATCH --time=06-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem=30gb								# Job memory request.

################################################################################
#                       02-Bismark.sh
#
#   This script runs the Bismark pipeline for whole-genome bisulfite sequencing
#   (WGBS) data processing. It includes:
#     - Genome preparation (indexing)
#     - Paired-end alignment
#     - Deduplication of BAM files
#     - Methylation extraction
#     - Global quality control and report generation
#
#   Input  : Raw paired-end reads (FASTQ)
#   Output : Aligned BAM files, methylation calls, and summary reports
#
################################################################################

####### MODULES
module load biotools

####### VARIABLES
WD="/home/nuezsal/Omics_integration/Methylome"
F="$WD/Scripts/02-Bismark/Functions.sh"

####### DIRECTORY
mkdir -p $WD/Results/02-Bismark
mkdir -p $WD/Results/02-Bismark/02-Alignment
mkdir -p $WD/Results/02-Bismark/03-Deduplication
mkdir -p $WD/Results/02-Bismark/04-Methylation_extractor
mkdir -p $WD/Results/02-Bismark/05-Global_reports
mkdir -p $WD/Results/02-Bismark/Logs

####### PIPELINE
#### BISMARK: GENOME PREPARATION (INDEXING)
echo -e "\n\n--------------------------------------------------"
echo -e "--------------- GENOME PREPARATION ---------------"
echo -e "--------------------------------------------------\n"

#srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Genome_preparation $WD/Results/02-Bismark

#### BISMARK: ALIGNMENT
echo -e "\n\n--------------------------------------------------"
echo -e "------------------- ALIGNMENT --------------------"
echo -e "--------------------------------------------------\n"


sample_list=$(ls $WD/Libraries/00-Raw_data)
for sample in $sample_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Alignment_paired_end $sample $WD/Libraries $WD/Results/02-Bismark &

done
wait

# Rename
sample_list=$(ls $WD/Libraries/00-Raw_data)
for sample in $sample_list; do
	mv $WD/Results/02-Bismark/02-Alignment/${sample}*_bismark_bt2_pe.bam $WD/Results/02-Bismark/02-Alignment/$sample.trimmed.bismark_bt2_pe.bam
	mv $WD/Results/02-Bismark/02-Alignment/${sample}*_bismark_bt2_PE_report.txt $WD/Results/02-Bismark/02-Alignment/$sample.trimmed.bismark_bt2_PE_report.txt 

done
wait

#### BISMARK: DEDUPLICATE
echo -e "\n\n--------------------------------------------------"
echo -e "----------------- DEDUPLICATION ------------------"
echo -e "--------------------------------------------------\n"


sample_list=$(ls $WD/Libraries/00-Raw_data)
for sample in $sample_list; do

	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Deduplication_paired_end $sample $WD/Results/02-Bismark &
done
wait

#### BISMARK: METHYLATION EXTRACTOR
echo -e "\n\n--------------------------------------------------"
echo -e "------------- METHYLATION EXTRACTOR --------------"
echo -e "--------------------------------------------------\n"

sample_list=$(ls $WD/Libraries/00-Raw_data)
for sample in $sample_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_Methylation_extractor_paired_end $sample $WD/Results/02-Bismark &
done
wait

#### BISMARK: GLOBAL REPORT
echo -e "\n\n--------------------------------------------------"
echo -e "----------------- GLOBAL REPORT ------------------"
echo -e "--------------------------------------------------\n"


sample_list=$(ls $WD/Libraries/00-Raw_data)
for sample in $sample_list; do
	srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive $F task_QC $sample $WD/Results/02-Bismark &
done
wait

multiqc $WD/Results/02-Bismark/05-Global_reports/ -o $WD/Results/02-Bismark/05-Global_reports/ >> $WD/Results/02-Bismark/Logs/MultiQC_stdout.log 2>&1

echo -e "\n\n*************** END BISMARK SCRIPT ***************\n\n"


