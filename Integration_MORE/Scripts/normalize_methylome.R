# Integration Omics with StategRa

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
suppressMessages(library("methylKit"))
suppressMessages(library(genomation))
suppressMessages(library(tibble))
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

files <- dir(path_table,pattern = "CX_report", full.names = TRUE)
metadata <- read.table(paste0(path_metadata,"/metadata_methylome.tsv"), sep = "\t", header = TRUE)
metadata$MORE <- gsub("_Sample", "", metadata$MORE)
# Define colors for each stress condition
condition_colors <- c(
  "cold" = "#87CEFA",         # Light blue
  "drought" = "#FFA07A",      # Salmon
  "monosporascus" = "#FF69B4",# Hot pink
  "control" = "#98FB98",      # Light green
  "shortday" = "#FFD700"      # Gold
)
metadata$Color <- condition_colors[metadata$Condition]

correspondance <- list("genes" = "CpG", "upstream" = "CHH")
Meth_mvalues <- data.frame()
for (type in c("genes", "upstream")){
  for(file in files){
    name <- basename(file)
    first_part <- strsplit(name, "\\.")[[1]][1]
    sample <- metadata[metadata$Sample == first_part, "MORE"]
    # read the files to a methylRawList object: myobj
    myobj=methRead(file,
                   sample.id=sample,
                   assembly="CMelo",
                   pipeline = 'bismarkCytosineReport',
                   treatment= 0,
                   context= correspondance[[type]],
                   mincov = 10
    )
    
    
    myobj_list <- new("methylRawList", list(myobj), treatment = 0) 
    
    regions <- read.table(paste0(path_annot,"/CMelon_DHL92_v4_",type,".bed"), sep = "\t")
    
    library(GenomicRanges)
    
    regions_GR <- GRanges(
      seqnames = regions$V1, 
      ranges = IRanges(start = regions$V2, end = regions$V3),
      strand = regions$V4,  # Añadir strand
      gene_id = regions$V5  # Añadir ID del gen como metadato
    )
    region_methyl <- regionCounts(myobj_list, regions_GR)
   
    data <- getData(region_methyl[[1]])
    # Combine chr, start, end, and strand to create the rowname column
    data$rowname <- paste(data$chr, data$start, data$end, ifelse(data$strand == "+", "F", "R"), sep = "_")
    
    # Calculate the 'sample' column as log2((numCs + 1) / (numTs + 1))
    data$sample <- log2((data$numCs + 1) / (data$numTs + 1))
    
    # Keep only the rowname and sample columns
    final_data <- unique(data[, c("rowname", "sample")])
    colnames(final_data) <- c("Position",sample)
    
    # Merge this final_data with Meth_mvalues based on the 'rownames'
    if (nrow(Meth_mvalues) == 0) {
      # If Meth_mvalues is empty, just assign final_data to it
      Meth_mvalues <- final_data
    } else {
      # Merge based on rownames, keeping all existing rows in Meth_mvalues
      Meth_mvalues <- merge(Meth_mvalues, final_data, by = "Position", all = FALSE)
    }
  }
  
  # Supongamos que 'df' es tu data frame y quieres que la columna 'ID' sea rownames
  Meth_mvalues_bp <- column_to_rownames(Meth_mvalues, var = "Position")
  
  boxplot(Meth_mvalues_bp,
          ylab = "M-value", 
          main = "Distribution of M-value in methylated regions",
          col = metadata$Color)
  png(paste0(path_graph_out, "/boxplot_filtered_counts_",type,".png"), width = 800, height = 600)
  dev.off()
  
  write.table(file = paste0(path_table_out,"/methylation_normalize_counts_",type,".tsv"),Meth_mvalues, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}


