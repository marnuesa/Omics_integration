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
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(scales))
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

# Create output paths
path_out_feature <- paste(path_out, '01-Feature', sep = '/')
path_out_global <- paste(path_out, '02-Global', sep = '/')
path_out_pie <- paste(path_out, '03-Pie', sep = '/')

# Create directories if they do not exist
dir.create(path_out_feature , recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_global, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_pie, recursive = TRUE, showWarnings = FALSE)

print("The arguments are correct, the analysis will start now...")

# Analysis of methylation profiles across various stress conditions and time points
times <- c("T1","T2","T3")
stresses <- c("C","D","MON","SD")

# Initialize final results table
final_table <- data.frame()
for (stress in stresses){
  for (time in times) {
    
        print(paste("The", stress, "-", time, "file is being analyzed..."))
    
        # Read input files for each methylation context
        CG_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CG_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CG_file) <- c("Chr", "Start", "End", "Methylation_type", "Proportion_diff", "ID", "Feature")
        CHG_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CHG_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CHG_file) <- c("Chr", "Start", "End", "Methylation_type","Proportion_diff","ID", "Feature")
        CHH_file <- read.table(paste0(path_in,"/",time,"-",stress,"_DMRs_Bins_CHH_genes.bed"),
                              sep = '\t', header = FALSE)
        colnames(CHH_file) <- c("Chr", "Start", "End", "Methylation_type", "Proportion_diff", "ID", "Feature")
        
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

final_table[is.na(final_table)] <- 0
write.table(final_table,file=paste0(path_out_global,"/Proportions_table.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)

print("Creating the specific feature graphs...")

# Data transformation for visualization
df_long <- final_table %>%
      pivot_longer(cols = starts_with("Proportion"), 
                   names_to = "Type", 
                   values_to = "Proportion") %>%
      mutate(Type = gsub("Proportion_", "", Type))

# Delete not necessary columns
df_long_filter <- df_long  %>%
  dplyr::select(-Count_CG,-Count_CHG,-Count_CHH) 

# Create two columns for proportion which depend on meth type
df_wide <- df_long_filter %>%
  pivot_wider(
    names_from = Meth_type,   
    values_from = Proportion, 
    names_prefix = "Proportion_" 
  )

# Generate bar plots for features with significant methylation changes divide by meth type
for(feature in unique(df_wide$Feature)){
  filtered_data <- df_long[df_long$Feature == feature ,]
  max_prop <- max(filtered_data$Proportion)
  if (max_prop > 0.005){
    # Filter
    filtered_wide <-  df_wide[df_wide$Feature == feature ,]
    filtered_wide$Proportion_Hypo <- -(filtered_wide$Proportion_Hypo)

    # Plotting
    plot <- ggplot(data = filtered_wide, mapping = aes(y = Stress, alpha = Time)) +
      geom_col(aes(x = Proportion_Hyper, fill = Stress), position = "dodge") +
      geom_col(aes(x = Proportion_Hypo, fill = Stress), position = "dodge") +
      facet_wrap(~ Type, ncol = 1, scales = "free_y") +
      scale_alpha_manual(values = c(0.4, 0.7, 1)) +
      scale_x_continuous(
        limits = c(-max_prop, max_prop),  # Agregar un pequeño espacio adicional
        breaks = c(-max_prop, 0, max_prop),  # Establecer los puntos de corte en los extremos y en el centro
        labels = function(x) {scales::number_format(accuracy = 0.001)(abs(x))}
      ) +
      annotation_custom(
        grob = grid::rectGrob(gp = grid::gpar(col = "white", fill = "white")), 
        xmin = -0.005, xmax = 0.005, ymin = -Inf, ymax = Inf)+
      labs(title = paste0("Proporciones de Metilación (Hyper vs Hyppo) in ", feature), x = "Proporción", y = "Estrés") +
      theme_minimal() +
      theme(legend.position = "bottom")

    # Save plot
    ggsave(paste0(path_out_feature,"/",feature,"_proportion.png"), 
           plot = plot, width = 20, height = 15,bg =" white") 
  }
}

print("Creating the global graphs...")

# Delete not necessary columns
table_new <- final_table[,-c(2,4,6,8)]

# Summarize the counts of different methylation contexts (CG, CHG, CHH) by stress level
df_sum <- table_new %>%
  group_by(Feature, Stress, Time) %>%
  summarise(
    Count_CG = sum(Count_CG),
    Count_CHG = sum(Count_CHG),
    Count_CHH = sum(Count_CHH),
    .groups = "drop"  # Eliminar el agrupamiento después de la operación
  )

# Calculate new proportions
df_proportion <- df_sum %>%
  group_by(Stress, Time) %>%
  mutate(
    total_CG = sum(Count_CG),
    total_CHG = sum(Count_CHG),
    total_CHH = sum(Count_CHH),
    Proportion_CG = Count_CG / total_CG,
    Proportion_CHG = Count_CHG / total_CHG,
    Proportion_CHH = Count_CHH / total_CHH
  ) %>%
  ungroup() %>%
  dplyr::select(Feature, Count_CG, Count_CHG, Count_CHH, Stress, Time, Proportion_CG, Proportion_CHG, Proportion_CHH)

# Create a long format of the data for plotting
final_long <- df_proportion %>%
  pivot_longer(cols = starts_with("Proportion_"), names_to = "Context", 
               values_to = "Proportion") %>%
  mutate(Context = gsub("Proportion_", "", Context)) %>%
  select(-contains("Count"))

# Generate adecuate labels
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
  # Per feature plot
  ggplot(plot_data, aes(x = Proportion, y = Time, fill = Feature)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Context , scales = "free_y", ncol = 1, strip.position = "left",) +
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

  # Save plot
  ggsave(plot = last_plot(), filename = paste0(path_out_global, "/Feature_analysis_",stress,".svg"), height = 10, width = 17, )
}

############################################## Context PIE CHARTS #####################################################
print("Creating the pie charts...")
# Convert the summarized data to a new long format for easier plotting pie chart
df_pie_long <- df_sum %>%
  pivot_longer(cols = starts_with("Count_"), names_to = "Context", values_to = "Count") %>%
  mutate(Context = gsub("Count_", "", Context))  # Remove "Count_" prefix from context names

# Compute the total counts per stress level
summary_counts <- df_pie_long %>%
  group_by(Stress) %>%
  summarise(total_counts = sum(Count))

# Merge total counts with the original data and compute percentages
df_pie_final <- df_pie_long %>%
  left_join(summary_counts, by = "Stress") %>%
  group_by(Stress) %>%
  mutate(percentage = Count / sum(Count))

# Create pie charts for each stress level
pie_charts <- list()
stress_levels <- unique(df_pie_final$Stress)

for (level in stress_levels) {
  # Filter data for the current stress level
  df_filtered <- df_pie_final %>% filter(Stress == level)
  
  # Ensure the Context variable is treated as a factor
  df_pie_final$Context <- factor(df_pie_final$Context)
  
  # Create a pie chart
  p <- ggplot(df_filtered, aes(x = "", y = percentage, fill = Context)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +  # Add white border to segments
    coord_polar("y", start = 0) +  # Convert to a pie chart
    scale_fill_manual(values = c("CG" = "#cb7eff", "CHG" = "#946eff", "CHH" = "#5d5dff")) +  # Custom colors
    labs(title = paste0(level, " - Total DMRs: ", unique(df_filtered$total_counts)), 
         fill = "Context") +
    theme_void() +  # Remove background and gridlines
    geom_text(aes(label = percent(percentage, accuracy = 0.1)),  # Format percentages
              position = position_stack(vjust = 0.5),
              color = "black",  # Black text for readability
              size = 3) +
    theme(legend.position = "bottom")  # Move legend to the bottom
  
  # Store the plot in a list
  pie_charts[[level]] <- p
}

# Define colors and labels for the legend
colors <- c("#cb7eff", "#946eff", "#5d5dff")  # Defined colors
contexts <- c("CG", "CHG", "CHH")  # Methylation contexts

# Create a legend as a graphical object
legend_manual <- legendGrob(
  labels = contexts, 
  pch = 15,  # Square symbols for legend
  gp = gpar(col = colors, fill = colors, fontsize = 12),
  ncol = length(contexts)  # Arrange in a single row (horizontal legend)
)

# Remove legends from individual pie charts
pie_charts_wo_legend <- lapply(pie_charts, function(g) g + theme(legend.position = "none"))

# Save the final image with multiple pie charts and a legend
png(paste0(path_out_pie,"/distribution_context.png"), 
    width = 10, height = 6, units = "in", res = 300)

# Arrange the pie charts in a grid and add the legend at the bottom
grid.arrange(do.call(arrangeGrob, c(pie_charts_wo_legend, ncol = 2, nrow = 2)), 
             legend_manual, 
             nrow = 2, heights = c(10, 1))

# Close the graphical device
dev.off()
