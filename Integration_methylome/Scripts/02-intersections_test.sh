#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=16               # Number of cpus
#SBATCH -t 0-01:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=30G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o intersections_%j.log # Standard output goes to this file

################################################################################
#
#               02-intersections_test.sh
#   
#    This program checks if there is significant overlap between the different
#    BED files.
#
################################################################################


path_out='/home/nuezsal/Omics_integration/Integration_methylome/Results/Intersections_test'
path_in='/home/nuezsal/Omics_integration/Integration_methylome/Results/Bed_files'

mkdir -p "$path_out"

bedtools intersect -wa -wb -a $path_in/CMelon_DHL92_v4_retrotransposons.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	 -filenames > $path_out/retrotransposons.txt

bedtools intersect -wa -wb -a $path_in/CMelon_DHL92_v4_downstream.bed \
	-b $path_in/CMelon_DHL92_v4_retrotransposons.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/genedowns.txt

bedtools intersect -wa -wb -a $path_in/CMelon_DHL92_v4_genes.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_retrotransposons.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/genes.txt

bedtools intersect -wa -wb -a $path_in/CMelon_DHL92_v4_upstream.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_retrotransposons.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/geneups.txt

bedtools intersect -wa -wb -a $path_in/microRNA_3_prime.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/CMelon_DHL92_v4_retrotransposons.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/micro3.txt

bedtools intersect -wa -wb -a $path_in/microRNA_5_prime.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed $path_in/microRNA_3_prime.bed \
	$path_in/CMelon_DHL92_v4_retrotransposons.bed  $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/micro5.txt

bedtools intersect -wa -wb -a $path_in/microRNA_precursors.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/CMelon_DHL92_v4_retrotransposons.bed $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-filenames > $path_out/micro.txt

bedtools intersect -wa -wb -a $path_in/CMelon_DHL92_v4_lncRNA.bed \
	-b $path_in/CMelon_DHL92_v4_downstream.bed $path_in/CMelon_DHL92_v4_genes.bed $path_in/CMelon_DHL92_v4_upstream.bed \
	$path_in/microRNA_3_prime.bed $path_in/microRNA_5_prime.bed $path_in/microRNA_precursors.bed $path_in/CMelon_DHL92_v4_retrotransposons.bed \
	-filenames > $path_out/lnc.txt
