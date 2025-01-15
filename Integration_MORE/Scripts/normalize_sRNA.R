################################################################################
##                                                                            
##  Normalization and filter sRNA count tables                 
##                                                                            
##  1. Metadata and Counts Processing
##
##  This program processes the metadata and counts data from an sRNA sequencing 
##  experiment. It normalizes the count data by DESeq2’s median of ratios, 
##  assigns stress-specific colors 
##  for visualization.
##
##  2. Filter by DE microRNA
##
##  The table with normalize counts is filtered by sequences which are annotated
##  as microRNA and have differential expression in some combination of time and 
##  stress. After that, all values of sequences that have the same microRNA annotation
##  are summed.
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
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library("argparse"))

#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'Normalize_sRNA.R',
                           description = '
    This program takes the tables of sRNA absolute counts and metadata:
    1. Metadata and Counts Processing
    
     This program processes the metadata and counts data from an sRNA sequencing
     experiment. It normalizes the count data by DESeq2’s median of ratios,
     assigns stress-specific colors
     for visualization.
    
     2. Filter by DE microRNA
    
     The table with normalize counts is filtered by sequences which are annotated
     as microRNA and have differential expression in some combination of time and
     stress. After that, all values of sequences that have the same microRNA annotation
     are summed.
    
     3. Data Visualization
    
     The script generates boxplots to visualize the distribution of raw and normalized
     counts, color-coded by stress conditions',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'sRNA project directory path.',
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
                        help = "path were DE microRNA are annotated")

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
path_table_out <- paste(path_out, '01-Tables', sep = '/')
path_graph_out <- paste(path_out, '02-Boxplots', sep = '/')

# Create directories if they do not exist
dir.create(path_table_out, recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_out, recursive = TRUE, showWarnings = FALSE)

print("######## Results files created ##########")

# Load count data and metadata
count_data <- read.table(paste0(path_table,"/fusion_abs-outer.csv"), 
                         sep = ",", header = TRUE, row.names = 1)
metadata <- read.table(paste0(path_metadata,"/metadata_sRNA.tsv"), 
                       sep = "\t", header = TRUE, stringsAsFactors = TRUE)

print("######## Count data loaded ##########")

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

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ Group)

# Pre-filtering.
keep <- rowSums(counts(dds) > 5) >= 10
dds<- dds[keep,]

print("######## Indepences filter done ##########")

# Raw count boxplot
raw_counts <- counts(dds)

png(paste0(path_graph_out, "/boxplot_raw_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((raw_counts + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Raw Counts",
        col = metadata$Color)
dev.off()

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

normalized_counts[is.na(normalized_counts)] <- 0  # Replace NAs with 0
normalized_counts[normalized_counts == 0] <- 1  # Replace zero values (to avoid log(0))

# Normalized count boxplot
png(paste0(path_graph_out, "/boxplot_normalized_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((normalized_counts + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Normalized Counts",
        col = metadata$Color)
dev.off()

print("######## Counts normalized ##########")

# Load and process annotated tables
annotated_files <- list.files(path_annot, full.names = TRUE)
annotated_results <- list()

for (i in seq_along(annotated_files)) {
  annotation <- read.table(annotated_files[i], header = TRUE, sep = ",")
  selected_columns <- annotation[, c(1, ncol(annotation))]
  annotated_results[[i]] <- unique(selected_columns)
}

combined_annotations <- do.call(rbind, annotated_results)
unique_annotations <- unique(combined_annotations)

# Filter normalized counts by annotations
filtered_counts <- normalized_counts[row.names(normalized_counts) %in% unique_annotations$seq, ]
filtered_counts_df <- filtered_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "seq")
merged_table <- merge(filtered_counts_df, unique_annotations, by = "seq")

# Handle duplicates
duplicates <- duplicated(merged_table$Row.names)
duplicate_rows <- merged_table[duplicates, ]
filtered_table <- merged_table[!merged_table$Row.names %in% duplicate_rows$Row.names, ]

# Summarize counts by annotation
final_table <- filtered_table[, -1] %>%
  group_by_at(ncol(.)) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var = "general_annot")

# Boxplot for filtered and summarized counts
png(paste0(path_graph_out, "/boxplot_filtered_counts.png"), width = 800, height = 600)
par(las = 2, mar = c(8, 5, 4, 2) + 0.1)
boxplot((final_table + 1), log = "y",
        ylab = "Log10(Normalized Counts + 1)", 
        main = "Distribution of sRNA Filtered Counts",
        col = metadata$Color)
dev.off()

# Save normalize count table
summed_table_final$miRNAs <- row.names(final_table)
write.table(file = paste0(path_table_out,"/sRNA_normalize_counts.tsv"),final_table,sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)
