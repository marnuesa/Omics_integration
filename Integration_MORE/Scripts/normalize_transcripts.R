################################################################################
##                                                                            
##  Normalization and filter transcripts count tables                 
##                                                                            
##  1. Metadata and Counts Processing
##
##  This program processes the metadata and counts data from a RNA-seq 
##  experiment. It normalizes the count data by DESeq2’s median of ratios, 
##  assigns stress-specific colors 
##  for visualization.
##
##  2. Filter by DE genes
##
##  The table with normalize counts is filtered by genes which are 
##  differential expression in some combination of time and stress. 
##
##  3. Data Visualization
##
##  The script generates boxplots to visualize the distribution of raw and normalized 
##  counts, color-coded by stress conditions.
##                                                                            
##  Author: Marta Núñez Salvador
##  Version: 1.0                                          
##  Date: 20/12/2024                                                          
##                                                                            
################################################################################

# Load necessary libraries
suppressMessages(library(DESeq2))
suppressMessages(library(tidyr))
suppressMessages(library(tximport))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library("argparse"))

#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'Normalize_transcripts.R',
                           description = '
    This program takes the tables of genes absolute counts and metadata:
    1. Metadata and Counts Processing
    
     This program processes the metadata and counts data from RNA-seq
     experiment. It normalizes the count data by DESeq2’s median of ratios,
     assigns stress-specific colors
     for visualization.
    
     2. Filter by DE microRNA
    
     The table with normalize counts is filtered by genes which are differential
     expression in some combination of time and stress.
    
     3. Data Visualization
    
     The script generates boxplots to visualize the distribution of raw and normalized
     counts, color-coded by stress conditions',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'RNAseq project directory path.',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Results path',
                        required = TRUE)
  required$add_argument('-m', '--metadata',
                        type = 'character',
                        help = 'path were metadata is save',
                        required = TRUE)
  required$add_argument('-a', '--annotations',
                        type = 'character',
                        help = "path were DE genes are save")
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('input', 'output', 'metadata','annotations')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  # Check if the input directory exists
  if (!dir.exists(args$input)) {
    stop('Error. The input directory does not exist.')
  }
  
  return(args)
}

##################################### MAIN #####################################

# Get programm arguments
args <- get_arguments()

# Save the the arguments in variables
path_table <- args$input
path_metadata <- args$metadata
path_out <- args$output
path_annot <- args$annotations

# Create output paths
path_table_out <- paste(path_out, '01-Tables',sep = '/')
path_graph_out <- paste(path_out, '02-Boxplots', sep = '/')

# Create directories if they do not exist
dir.create(path_table_out, recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_out, recursive = TRUE, showWarnings = FALSE)

print("######## Results files created ##########")
       
# Load count data and metadata
DE_genes <- read.table(paste0(path_annot,"/DE_genes.txt"))
metadata <- read.table(paste0(path_metadata,"/metadata_transcripts.tsv"), 
                       sep = "\t", header = TRUE, stringsAsFactors = TRUE)
# Clean and prepare metadata
metadata$SampleID <- gsub("_Sample", "", metadata$MORE)

# Define colors for each stress condition
condition_colors <- c(
  "cold" = "#87CEFA",         # Light blue
  "drought" = "#FFA07A",      # Salmon
  "monosporascus" = "#FF69B4",# Hot pink
  "control" = "#98FB98",      # Light green
  "shortday" = "#FFD700"      # Gold
)
metadata$Color <- condition_colors[metadata$Condition]

files <- file.path(path_table, metadata$Sample, "quant.sf")
names(files) <- row.names(metadata)

# The tx2gene file was created using the transcriptome file headers.
tx2gene <- read.table(paste0(path_metadata,"/tx2gene.txt"), sep = "\t", header = FALSE, stringsAsFactors = TRUE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no") 

# Let's construct a DESeqDataSet from the txi `object` and sample information in `sampletable`
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = metadata,
                                   design = ~Group)
# Pre-filtering.
keep <- rowSums(counts(ddsTxi) > 5) >= 10
ddsTxi<- ddsTxi[keep,]

# Raw count boxplot
raw_counts <- counts(ddsTxi)

png(paste0(path_graph_out, "/boxplot_raw_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((raw_counts + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Raw Counts",
        col = metadata$Color)
dev.off()

# Normalize counts
dds <- estimateSizeFactors(ddsTxi)
normalized_counts <- counts(dds, normalized = TRUE)

print("######## Counts normalized ##########")
       
# Normalized count boxplot
png(paste0(path_graph_out, "/boxplot_normalized_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((normalized_counts + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Normalized Counts",
        col = metadata$Color)
dev.off()

normalized_counts_filter <- normalized_counts[row.names(normalized_counts) %in% DE_genes$V1,]

# Boxplot for filtered and summarized counts
png(paste0(path_graph_out, "/boxplot_filtered_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((final_table + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Filtered Counts",
        col = metadata$Color)
dev.off()

# Save normalize count table
normalized_counts_filter$genes <- row.names(normalized_counts_filter)
write.table(file = paste0(path_table_out,"/Genes_normalize_counts.tsv"),normalized_counts_filter,sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)
