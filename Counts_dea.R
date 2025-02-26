# ==========================================================
#                    Counts_dea.R
# ==========================================================
# Description:
# ----------------------------------------------------------
# This script processes omics data from three different sources:
# 1. Methylome (DMRs)
# 2. MicroRNA
# 3. RNA-seq
# 
# It calculates the number of differentially expressed or regulated elements
# for each omic type and organizes the data across various time points and 
# stress conditions.
# 
# 
# Author: Marta Núñez Salvador
# Date: 22/09/2024
# Version: 1.1
################################################################################

suppressMessages(library(argparse))
suppressMessages(library(dplyr))

################################## FUNCTIONS ###################################
# Function to parse command-line arguments
get_arguments <- function() {
  parser <- ArgumentParser(prog = 'Conteos_dea.R',
                           description = 'Processes omics data and generates summary tables.')
  
  required <- parser$add_argument_group('required arguments')
  
  required$add_argument('-m', '--methylome_path', type = 'character',
                        help = 'Path to the methylome data directory.', required = TRUE)
  required$add_argument('-mi', '--microrna_path', type = 'character',
                        help = 'Path to the microRNA data directory.', required = TRUE)
  required$add_argument('-r', '--rnaseq_path', type = 'character',
                        help = 'Path to the RNA-seq data directory.', required = TRUE)
  required$add_argument('-o', '--output_path', type = 'character',
                        help = 'Path to the output directory.', required = TRUE)
  
  args <- parser$parse_args()
  return(args)
}

################################# MAIN #########################################
# Get program arguments
args <- get_arguments()

# Save output arguments in variables
output_path <- args$output_path

# Create output directory if it does not exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

#################### Methylome ####################
# Define the path for methylome data
path_methylome <- args$methylome_path

# List all the methylome TSV files
files_methylome <- list.files(path_methylome, pattern = ".tsv", full.names = TRUE)

# Initialize an empty dataframe for methylome
df_methylome <- data.frame()
for(file in files_methylome) {
  # Extract the filename from the full file path
  file_name <- basename(file)
  
  # Split the filename using '.' and '_' as delimiters, and extract the stress condition
  stress <- strsplit(file_name, "[._]")[[1]][3]
  
  # Read the TSV file into a data frame
  table <- read.table(file, header = TRUE, sep = "\t")
  
  # Add a new column called 'Stress' with the extracted value
  table$Stress <- stress
  
  # Combine the current file's data with the main methylome data frame
  df_methylome <- rbind(df_methylome, table)
}

# Convert the Time column to a factor if it's not already
df_methylome$Time <- as.factor(df_methylome$Time)

#################### MicroRNA ####################
# Define the path for microRNA data
path_microrna <- args$microrna_path

# List all microRNA files
files_microrna <- list.files(path_microrna, full.names = TRUE)

# Initialize an empty list for microRNA data
df_microrna <- list()
for (file in files_microrna){
  # Extract the filename from the full file path
  file_name <- basename(file)
  
  # Split the filename using '.' and '_' as delimiters, and extract stress and time
  stress <- strsplit(file_name, "[._]")[[1]][1]
  time <- strsplit(file_name,"[._]")[[1]][2]
  
  # Read the CSV file into a data frame
  table <- read.table(file, sep = ",", header = TRUE)
  
  # Add a column for differential expression (up or down-regulated)
  table$DE_type <- ifelse(table$Shrunkenlog2FoldChange < 0, "Infra-expresado", "Sobre-expresado")
  
  # Select only relevant columns and add Stress and Time columns
  table <- table[, c("seq", "DE_type")]
  table$Stress <- stress
  table$Time <- time
  
  # Combine the current file's data with the main microRNA data frame
  df_microrna <- rbind(df_microrna, table)
}

# Group by Time, Stress, and DE_type to count occurrences
df_microrna_counts <- df_microrna %>%
  group_by(Time, Stress, DE_type) %>%
  summarise(Count = n(), .groups = 'drop')

#################### RNA-seq ####################
# Define the path for RNA-seq data
path_rnaseq <- args$rnaseq_path

# List all RNA-seq files
files_rnaseq <- list.files(path_rnaseq, full.names = TRUE)

# Initialize an empty list for RNA-seq data
df_rnaseq <- list()
for (file in files_rnaseq){
  # Extract the filename from the full file path
  file_name <- basename(file)
  
  # Split the filename using '.' and '_' as delimiters, and extract stress and time
  stress <- strsplit(file_name, "[._]")[[1]][1]
  time <- strsplit(file_name,"[._]")[[1]][2]
  
  # Read the CSV file into a data frame
  table <- read.table(file, sep = ",", header = TRUE)
  if (nrow(table) != 0) {
    # Add a column for differential expression (up or down-regulated)
    table$DE_type <- ifelse(table$Shrunkenlog2FoldChange < 0, "Infra-expresado", "Sobre-expresado")
    
    # Select only relevant columns and add Stress and Time columns
    table <- table[, c("seq", "DE_type")]
    table$Stress <- stress
    table$Time <- time
    
    # Combine the current file's data with the main RNA-seq data frame
    df_rnaseq <- rbind(df_rnaseq, table)
  }
}

# Group by Time, Stress, and DE_type to count occurrences
df_rnaseq_counts <- df_rnaseq %>%
  group_by(Time, Stress, DE_type) %>%
  summarise(Count = n(), .groups = 'drop')


write.table(df_methylome,file = paste0(output_path,"/summary_methylome.tsv"),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(df_microrna_counts,file = paste0(output_path,"/summary_microRNA.tsv"),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(df_rnaseq_counts,file = paste0(output_path,"/summary_transcripts.tsv"),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
