# ==========================================================
#                      Plots.R
# ==========================================================
# 
# Description:
# ----------------------------------------------------------
# This script processes methylation data across various stress conditions
# (Cold, Drought, Monosporascus, Short Day) and time points to analyze 
# differentially methylated regions (DMRs) in different methylation 
# contexts (CG, CHG, CHH). The script generates visualizations of the 
# distribution of hypo- and hyper-methylated regions for each condition.
#
# Author: Marta Nuñez Salvador
# Date: 29/08/24
# Version: 1.0
#
#
# Notes:
# ----------------------------------------------------------
# - This script reads DMR files located in the specified directory and 
#   categorizes them by stress conditions using pattern matching.
# - The input files should be tab-delimited with columns such as Chr, 
#   Start, End, Methylation_type, ID, and Feature.
# - The methylation type is standardized from "gain" to "HyperMethylated" 
#   and "loss" to "HypoMethylated".
# - The script calculates the number of hypo- and hyper-methylated regions 
#   for each context (CG, CHG, CHH) and time point.
# - Visualization is done through bar plots showing the counts of 
#   methylated regions across different contexts and stress conditions.
#
# Output:
# ----------------------------------------------------------
# - Bar plots illustrating the number of hypo- and hyper-methylated 
#   regions for each context (CG, CHG, CHH) and stress condition.
# - The plots are saved as SVG files in the specified output directory.
# - A TSV file for each stress condition, containing the counts of hypo- 
#   and hyper-methylated regions.
# ==========================================================

rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(argparse))
suppressMessages(library(gridExtra))
################################## FUNCTIONS ###################################
# Get the command line arguments
# This function parse the command line arguments entered into the program.
#
# @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'DMRs_plot.R',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'Input directory path.',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Output directory path',
                        required = TRUE)
  
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('input', 'output')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  return(args)
}

# Function to create sublists based on a pattern
create_sublist <- function(rutas, patron) {
  sublist <- rutas[grep(patron, rutas)]
  
  if (length(sublist) == 0) {
    print(paste("No matches found for the pattern:", patron))
  }
  return(sublist)
}

################################# MAIN #########################################
# Get programm arguments
args <- get_arguments()

# Save the the arguments in variables
directory <- args$input
path_out <- args$output

print("The arguments are correct, the analysis will start now...")

# Get the list of files in the directory
files <- list.files(path = directory, full.names = TRUE)

# Create sublists for each category
files_C <- create_sublist(files, "T[0-9]-C.*duplicates_removed\\.tsv")
files_SD <- create_sublist(files, "T[0-9]-SD.*duplicates_removed\\.tsv")
files_MON <- create_sublist(files, "T[0-9]-MON.*duplicates_removed\\.tsv")
files_D <- create_sublist(files, "T[0-9]-D.*duplicates_removed\\.tsv")

stresses_list <- list(files_C, files_D, files_MON, files_SD)

# Add stress type to each sublist
stress_types <- c("Cold", "Drought", "Monosporascus", "Short Day")
for (i in seq_along(stresses_list)) {
  estres <- stress_types[i]
  names(stresses_list)[i] <- estres
}

plots <- list()
for (estres in names(stresses_list)) {
  files_list <- stresses_list[[estres]]
  all_counts <- data.frame()
  
  for(file in files_list){
    
    print(paste("The", file,"is being analyzed..."))
    
    # Read file
    stress_table <- read.table(file, sep = '\t', header = TRUE)
    
    # Extract hypo and hyper variables
    stress_table$regionType[stress_table$regionType == "gain"] <- "HyperMethylated"
    stress_table$regionType[stress_table$regionType == "loss"] <- "HypoMethylated"
    counts <- table(stress_table$regionType)
    
    # Generate columns
    name_parts <- unlist(strsplit(basename(file), "-|_"))
    time <- name_parts[1]
    context <- name_parts[5]
    hypo <- as.numeric(counts["HypoMethylated"])
    hyper <- as.numeric(counts["HyperMethylated"])
    
    # Generate rows
    new_row_hyper <- data.frame(Context = context,
                                   Time = time,
                                   Value = hyper,
                                   MethylatedType = "HyperMethylated")
    new_row_hypo <- data.frame(Context = context,
                                  Time = time,
                                  Value = hypo,
                                  MethylatedType = "HypoMethylated")
    
    # Add to the dataframe
    all_counts <- rbind(all_counts, new_row_hyper)
    all_counts <- rbind(all_counts, new_row_hypo)
  }
  
  print("Creating plots...")
	
  # Saving table
  write.table(all_counts, file = paste0(path_out, "/all_counts_", estres, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)  
  
  # Create the first plot for the mean and standard deviation of CG
  plot1 <- ggplot(all_counts[all_counts$Context == "CG", ], aes(x = Time, y = Value, fill = MethylatedType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "CpG", y = "Number DMRs", x = "Time") +
    scale_fill_manual(values = c("HyperMethylated" = "#00aaff", "HypoMethylated" = "#f34336")) +
    ylim(0, 130000) +
    theme_minimal()
  
  # Create the second plot for the mean and standard deviation of CHG
  plot2 <- ggplot(all_counts[all_counts$Context == "CHG", ], aes(x = Time, y = Value, fill = MethylatedType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "CHG", y = "Number DMRs", x = "Time") +
    scale_fill_manual(values = c("HyperMethylated" = "#00aaff", "HypoMethylated" = "#f34336")) +
    ylim(0, 130000) +
    theme_minimal()
  
  # Create the third plot for the mean and standard deviation of CHH
  plot3 <- ggplot(all_counts[all_counts$Context == "CHH", ], aes(x = Time, y = Value, fill = MethylatedType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "CHH", y = "Number DMRs", x = "Time") +
    scale_fill_manual(values = c("HyperMethylated" = "#00aaff", "HypoMethylated" = "#f34336")) +
    ylim(0, 130000) +
    theme_minimal()
  
  # Combine the three plots into one visualization
  plot <- grid.arrange(plot1, plot2, plot3, ncol = 3)
  plots[[estres]] <- plot
  
  print("Saving plots...")
  
  # Save plot
  ggsave(paste0(path_out,"/plot_", estres, ".svg"), plot, width = 20, height = 15, dpi = 600)

}
