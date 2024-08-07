#!/bin/bash
#SBATCH --ntasks=10               # Number of task
#SBATCH --cpus-per-task=16               # Number of cpus
#SBATCH -t 1-00:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=5G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o Trimming_%j.log # Standard output goes to this file

################################################################################
#
#               01-Trimming.sh
#   This function trims the sequences from a
#   specific library using the TrimGalore program:
#   1- Adapter trimming
#   2- poliA trimming
#
#   It is necessary to do it this way because otherwise the polyA is not
#   trimmed and the libraries are not cleaned well.
#
#   Arguments:
#       Fastq.gz file path
#       Output directory path
#
################################################################################
# MODULE
module load anaconda

# PATHs
path_in='/home/nuezsal/Omics_integration/Methylome/Libraries/00-Raw_data'
path_out_adapter='/home/nuezsal/Omics_integration/Methylome/Libraries/01-Clean_data_adapter'
path_out='/home/nuezsal/Omics_integration/Methylome/Libraries/01-Clean_data'

# Create output fields
if [ ! -d "$path_out_adapter" ]; then
    echo "El directorio de salida no existe. Creando directorio..."
    mkdir "$path_out_adapter"
    echo "Directorio creado correctamente."
else
    echo "El directorio de salida ya existe."
fi

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "El directorio de salida no existe. Creando directorio..."
    mkdir "$path_out"
    echo "Directorio creado correctamente."
else
    echo "El directorio de salida ya existe."
fi

# List species directories
sample_list=$( ls ${path_in})

## ADAPTER TRIMMING
# Iterate files list
#for sample in $sample_list
#do
    # Trimming
#    echo "Trimming $sample sample"
#    srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive trim_galore --paired \
#        --phred33 \
#        --quality 20 \
#        --length 20 \
#        --clip_R1 10 \
#        --clip_R2 10 \
#        --adapter GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#        --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
#        -o ${path_out_adapter} --basename $sample --cores 4 \
#        ${path_in}/${sample}/${sample}_1.fq.gz ${path_in}/${sample}/${sample}_2.fq.gz &

#done
#wait

## POLI A TRIMMING
# Iterate files list
for sample in $sample_list
do
    # Trimming
    echo "Trimming $sample sample"
    srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive trim_galore --paired \
        --phred33 \
        --quality 20 \
        --length 20 \
        --clip_R1 10 \
        --clip_R2 10 \
        -a A{10} -a2 A{10} \
        -o ${path_out} --basename $sample --cores 4 \
        ${path_out_adapter}/${sample}_val_1.fq.gz ${path_out_adapter}/${sample}_val_2.fq.gz &

done
wait
exit 0

# rm -r $path_out_adapter
