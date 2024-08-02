#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=6              # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem-per-cpu=20G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o 06-Creating_STARindex_%j.log # Standard output goes to this file

#VARIABLES
genoma_gff="/home/nuezsal/omics_integration_G3/RNAseq/Additional_data/Melon_PacBio_v4_liftover.gff"
genoma_fa="/home/nuezsal/omics_integration_G3/RNAseq/Additional_data/Melon_v4.0_PacBio.fasta"
index_dir="/home/nuezsal/omics_integration_G3/RNAseq/Results/02-Salmon/index_STAR"


#Create an index

srun STAR --runThreadN 6 \
	--runMode genomeGenerate \
	--genomeDir ${index_dir} \
	--genomeFastaFiles ${genoma_fa} \
	--sjdbGTFfile  ${genoma_gff} --sjdbGTFtagExonParentTranscript Parent \
	--sjdbOverhang 100
