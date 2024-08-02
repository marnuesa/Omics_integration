#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=10                # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=25G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o salmon_prep_%j.log  # Standard output goes to this file

################################################################################
#
#               01-Salmon_index.sh
#   This program create a Salmon index using the Melon genome v4 of Cucurbit 
#   genomic web v2.0
#
################################################################################

# MODULES
module load anaconda

# PATHs
path_genome=/home/nuezsal/Omics_integration/RNA_seq/Additional_info

# Use the entire genome as decoy sequence
grep "^>" "$path_genome/DHL92_genome_v4.fa" | cut -d " " -f 1 | sed 's/>//g' > "$path_genome/decoys.txt"

# Along with the list of decoys salmon also needs the concatenated transcriptome and 
# genome reference file for index
cat "$path_genome/DHL92_mRNA_v4.fa" "$path_genome/DHL92_genome_v4.fa" > "$path_genome/gentrome.fa"

# Create the salmon index
srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive salmon index -t "$path_genome/gentrome.fa" -i "$path_genome/Melon_salmon_ind" -d "$path_genome/decoys.txt"
