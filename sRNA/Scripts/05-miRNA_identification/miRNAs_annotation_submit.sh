#!/bin/bash

#******************************************************************************
#  
#   miRNAs_annotation_submit.sh
#
#   The program is designed to identify miRNAs from tables of differentially
#   expressed sequences for multiple species using the miRNAS_annotation.sh
#   programm.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 22/12/2023
#   Version: 2.0 
#
#******************************************************************************

# Input paths
path_in=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Results/04-deseq_sig
path_mirbase=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Additional_info/02-Mod_databases/miRBase
path_PmiREN=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Additional_info/02-Mod_databases/PmiREN
path_sRNAanno=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Additional_info/02-Mod_databases/sRNAanno
path_ids_table=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Additional_info/species_id.csv
mismatches=1

# Ouput paths
path_out=/home/marnuesa/Documentos/Omics_integration/sRNA_reverse/Results/05-Identification_miRNAs

# Threads
num_threads=12

# Execution 
bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mismatches $mismatches \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads 

exit 0
