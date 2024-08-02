#!/bin/bash
#SBATCH --ntasks=1               # Number of task
#SBATCH --cpus-per-task=16               # Number of cpus
#SBATCH -t 0-03:00:00       # Runtime in minutes.
#SBATCH --qos short         # The QoS to submit the job.
#SBATCH --mem=10G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o post-fastqc_%j.log # Standard output goes to this file

################################################################################
#                       00-PreFastQC.sh
#
#   This script execute fastqc to do a preliminar quality control with
#   Raw data
#
################################################################################

#MODULE
module load biotools

#PATHs
path_in='/home/nuezsal/Omics_integration/Methylome/Libraries/01-Clean_data'
path_out='/home/nuezsal/Omics_integration/Methylome/Results/01-PostFastQC'

# Create output fields
if [ ! -d "$path_out" ]; then
    echo "El directorio de salida no existe. Creando directorio..."
    mkdir "$path_out"
    echo "Directorio creado correctamente."
else
    echo "El directorio de salida ya existe."
fi

# Execute
srun -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive fastqc ${path_in}/*fq.gz -t 16 -o ${path_out}
