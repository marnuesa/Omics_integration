# ==========================================================
#                      04-DMRs_feature_plot.R
# ==========================================================
# 
# Description:
# ----------------------------------------------------------
# This script analyzes methylation profiles across various stress
# conditions and time points, focusing on different types of 
# methylation (CG, CHG, CHH) and their associated features.
# The aim is to identify patterns in methylation across different 
# contexts and stress conditions, and visualize the distribution 
# of differentially methylated regions (DMRs) and their features.
#
# Author: Marta Nuñez Salvador
# Date: 27/08/24
# Version: 1.0
#
#
# Notes:
# ----------------------------------------------------------
# - This script requires input files containing DMRs (Differentially 
#   Methylated Regions) for different methylation contexts (CG, CHG, CHH)
#   and conditions (time points and stress conditions).
# - Input files must be formatted correctly with columns: Chr, Start, End,
#   Methylation_type, ID, and Feature.
# - The methylation type is modified from "loss" to "Hypo" and other types 
#   to "Hyper".
# - Proportions of features for each methylation type are calculated and 
#   visualized.
# - Missing values are replaced with zeros before analysis.
#
# Output:
# ----------------------------------------------------------
# - Plots showing the distribution of DMRs across different contexts 
#   (CG, CHG, CHH) and stress conditions, with separate plots for each 
#   stress condition. Each plot displays the proportion of features in 
#   different contexts and time points, stacked by methylation type.
# - SVG files saved in the specified output directory.
# ==========================================================

rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(argparse))

################################## FUNCTIONS ###################################

# Define a function to modify the methylation type
modify_methylation_type <- function(df) {
  df$Methylation_type <- ifelse(df$Methylation_type == "loss", "Hypo", "Hyper")
  return(df)
}

# Define a function to calculate the proportion of each feature
calculate_proportion <- function(df) {
  prop_table <- as.data.frame(table(df$Feature, df$Methylation_type))
  colnames(prop_table) <- c("Feature", "Meth_type", "Count")
  prop_table <- prop_table %>%
    group_by(Meth_type) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  return(prop_table)
}

# Get the command line arguments
# This function parse the command line arguments entered into the program.
#
# @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = '04-DMRs_feature_plot.R',
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

################################# MAIN #########################################
# Get programm arguments
args <- get_arguments()

# Save the the arguments in variables
path_in <- args$input
path_out <- args$output

print("The arguments are correct, the analysis will start now...")

# Analysis of methylation profiles across various stress conditions and time points
times <- c("T1","T2","T3")
stresses <- c("C","D","MON","SD")
final_table <- data.frame()
for (stress in stresses){
  for (time in times) {
    
        print(paste("The", stress, "-", time, "file is being analyzed..."))
    
        CG_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CG_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CG_file) <- c("Chr", "Start", "End", "Methylation_type", "ID", "Feature")
        CHG_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CHG_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CHG_file) <- c("Chr", "Start", "End", "Methylation_type", "ID", "Feature")
        CHH_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CHH_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CHH_file) <- c("Chr", "Start", "End", "Methylation_type", "ID", "Feature")
        
        # Apply the function to modify the methylation type
        CG_file <- modify_methylation_type(CG_file)
        CHG_file <- modify_methylation_type(CHG_file)
        CHH_file <- modify_methylation_type(CHH_file)
        
        print("Calculating proportions...")
        
        # Calculate proportions for each data frame
        CG_proportions <- calculate_proportion(CG_file)
        CHG_proportions <- calculate_proportion(CHG_file)
        CHH_proportions <- calculate_proportion(CHH_file)
        
        print("Merge tables...")
        
        # Merge the proportion tables
        merged_proportions <- merge(CG_proportions, CHG_proportions, by = c("Feature", "Meth_type"), 
                                    all = TRUE, suffixes = c("_CG", "_CHG"))
        merged_proportions <- merge(merged_proportions, CHH_proportions, by = c("Feature", "Meth_type"), all = TRUE)
        colnames(merged_proportions)[colnames(merged_proportions) == "Proportion"] <- "Proportion_CHH"
        colnames(merged_proportions)[colnames(merged_proportions) == "Count"] <- "Count_CHH"
        
        # Add stress and time columns
        merged_proportions$Stress <- stress
        merged_proportions$Time <- time
        
        # Combine with the final table
        final_table <- rbind(final_table, merged_proportions)
  }
}

write.table(final_table,file=paste0(path_out,"/Proportions_table.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)

print("Creating the graphs...")

# Replace NA per 0 
final_table <- replace(final_table, is.na(final_table), 0)

# Create a long format of the data for plotting
final_long <- final_table %>%
  pivot_longer(cols = starts_with("Proportion_"), names_to = "Context", 
               values_to = "Proportion") %>%
  mutate(Context = gsub("Proportion_", "", Context)) %>%
  select(-contains("Count"))

final_long$Feature <- factor(final_long$Feature, levels = c("upstream","genes",
                                                            "5prime", "precursors",
                                                            "5primelncRNA", "lncRNA",
                                                            "downstream","3prime",
                                                            "3primelncRNA","retrotransposons",
                                                            "unknown_region"))
final_long$Time <- factor(final_long$Time, levels=c("T3","T2","T1"))
stress_correspondence <- list(C= "Cold", D="Drought",MON="Monosporascus",SD="Short Day")

# Plot
for (stress in unique(final_long$Stress)) {
  
  print(paste("Creating", stress, "plot"))
  plot_data <- final_long %>% filter(Stress == stress)

  ggplot(plot_data, aes(x = Proportion, y = Time, fill = Feature)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Context + Meth_type , scales = "free_y", ncol = 1, strip.position = "left",) +
    labs(title = paste("Distribution of DMRs in", stress_correspondence[[stress]]),
         x = "" ,
         y = "",
         fill = "Region") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 14),
          plot.title = element_text(family = "Helvetica", face = "bold", size = (17), hjust = 0.5)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_reverse() +
    scale_fill_manual(values = c("upstream" = "midnightblue","genes" = "mediumslateblue",
                                 "5prime"= "#0000CD", "precursors"="#00F5FF",
                                 "5primelncRNA"= "#87CEFA", "lncRNA"="#21D8AE",
                                 "downstream"="#98FB98","3prime"="orangered1", 
                                 "3primelncRNA"="#FFB90F", "retrotransposons"="#FFD39B",
                                 "unknown_region"="#FFF68F"), 
                      labels = c("upstream" = "Gene upstream","genes" = "Gene",
                                 "5prime"= "MicroRNA upstream", "precursors"="MicroRNA", 
                                 "5primelncRNA"= "lncRNA upstream","lncRNA"="lncRNA",
                                 "downstream"="Gene downstream","3prime"="MicroRNA downstream", 
                                 "3primelncRNA"="lncRNA_downstream","retrotransposons"="Retrotransposon",
                                 "unknown_region"="Unknown"))
  
  ggsave(plot = last_plot(), filename = paste0(path_out, "/Feature_analysis_",stress,".svg"), height = 10, width = 17, )
}

##################################################################################################################
