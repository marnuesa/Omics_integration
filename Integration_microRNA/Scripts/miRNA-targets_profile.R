# ==========================================================
#                          miRNA-targets_profile.R
# ==========================================================
# 
# Description:
# ----------------------------------------------------------
# This script analize the expression profile of microRNA
# which are differentially expressed at almost 3/4 stress conditions
# and their respetive targets. 
# The aim is identify a core of microRNA in respons to multiple
# stresses and look how their targets act.
#
# Author: Marta Núñez
# Data: 18/06/24
# Versión: 1.0
#
#
# Notes:
# ----------------------------------------------------------
# - This script needs the raw files of transcripts and microRNA
# DEA and the significant sequences of microRNA
# - microRNA has to be differentially significant in at least one 
# combination of stress/time to be considered.
# - microRNA hasn't to be eliminated for low counts in any time
# - Targets need to be differentially significant at list in one
# combination of stress/time and not be eliminated for low counts in any time.

# Output:
# ----------------------------------------------------------
# - Analysis 1: Diversity of microRNA sequences between stresses
# - Analysis 2: miRNA profile with high Basemean represents their miRNA family
# and is graphicated with each one of their targets profile. It generation a correlation
# plot, a heatmap and the edge and node tables to create a network.
# - Analysis 3: Select only the microRNAs that have a different profile that
# their family and and graphicated with each one of their targets profile.

rm(list = ls())

################################# Libraries ####################################
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))
suppressMessages(library("argparse"))
suppressMessages(library(svglite))
suppressMessages(library("cowplot"))
suppressMessages(library("plotly"))
suppressMessages(library(htmlwidgets))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library(dplyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.CMelo.eg.db))

################################## FUNCTIONS ###################################


#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'miRNA-targets_profile.R',
                           description = '
    This program takes the tables of absolute counts and',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-m', '--microrna',
                        type = 'character',
                        help = 'microRNA DEA directory path.',
                        required = TRUE)
  required$add_argument('-t', '--transcripts',
                        type = 'character',
                        help = 'transcripts raw DEA directory path',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'output path',
                        required = TRUE)
  required$add_argument('-ai', '--additional',
                      type = 'character',
                      help = 'Additional info path',
                      required = TRUE)
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('microrna','transcripts', 'output', 'additional')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  return(args)  
}
##################################### MAIN #####################################

# Get programm arguments
args <- get_arguments()

# Save the the arguments in variables
microRNA_path <- args$microrna
transcripts_raw_path <- args$transcripts
path_out <- args$output
path_ai <- args$additional

################################### PATHs ######################################

microRNA_raw_path <- paste0(microRNA_path, '/Group_miRNAs_raw/cume/Omics_project/01-DEA_results_annot')
microRNA_sig_path <- paste0(microRNA_path, '/Group_miRNAs_sig/cume/Omics_project/01-DEA_results_annot')

# Create output paths
path_out_analysis1 <- paste(path_out, '01-Analysis1', sep = '/')
path_out_analysis2 <- paste(path_out, '02-Analysis2', sep = '/')
path_out_analysis3 <- paste(path_out, '03-Analysis3', sep = '/')
path_out_logs <- paste(path_out, '04-Logs', sep = '/')
path_out_network <- paste(path_out, 'Network_files', sep = '/')

# Create directories if they do not exist
dir.create(path_out_analysis1 , recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_analysis2, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_analysis3, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_network, recursive = TRUE, showWarnings = FALSE)

################################## Read tables #################################

targets <- read.table(paste0(path_ai,"/Targets_orthologues.txt"),sep='\t', header = TRUE)
annotation_table <- read.table(paste0(path_ai,"/DHL92_gene_description_v4.txt"), sep = '\t',header=FALSE, quote = "")
colnames(annotation_table) <- c("gene", "description")

# Load significative sequences dataframe
microRNA_sig_files <- list.files(microRNA_sig_path, full.names = TRUE)
microRNA_raw_files <- list.files(microRNA_raw_path, full.names = TRUE)
transcripts_raw_files <- list.files(transcripts_raw_path, full.names = TRUE)

## Extract file names
sig_names <- sapply(strsplit(basename(microRNA_sig_files), "_"), `[`, 1)
raw_names <- sapply(strsplit(basename(microRNA_raw_files), "_"), `[`, 1)

## Create a list of directories organized in stresses
microRNA_sig <- split(microRNA_sig_files, sig_names)
microRNA_raw <- split(microRNA_raw_files, raw_names)
transcripts_raw <- split(transcripts_raw_files, raw_names)

# Load microRNA significant dataframes
dataframes_micro_sig <- list()
for(stress in names(microRNA_sig)){
  tablas <- lapply(microRNA_sig[[stress]],read_csv)
  dataframes_micro_sig[[stress]] <- tablas    
}

# Load microRNA and transcripts raw dataframes
dataframes_micro_raw <- list()
dataframes_trans <- list()
for(stress in names(microRNA_raw)){
  tablas_micro <- lapply(microRNA_raw[[stress]],read_csv)
  tablas_trans <- lapply(transcripts_raw[[stress]],read_csv)
  dataframes_micro_raw[[stress]] <- tablas_micro
  dataframes_trans[[stress]] <- tablas_trans
}

############################ Analysis 1  ########################################

# Extract microRNA sequences for stress
miRNAs <- list()
for(stress in names(dataframes_micro_sig)){
  micro_sig_stress <- dataframes_micro_sig[[stress]]
  micro_sig_stress_table <- do.call(rbind, micro_sig_stress)
  miRNAs_annot <- unique(micro_sig_stress_table$general_annot)
  
  for (miRNA in miRNAs_annot){
    miRNAs_seq <- micro_sig_stress_table[micro_sig_stress_table$general_annot==miRNA,]
    miRNAs_seq <- unique(miRNAs_seq$seq)
    miRNAs[[miRNA]][[stress]] <- miRNAs_seq
  }
}

# Datraframe to save microRNA
micro_summary_plot <- data.frame()
micro_summary_table <- data.frame()
patron = 0
total = 0
# Iterate microRNAs
for(miRNA in names(miRNAs)){
  miRNAs_list <- miRNAs[[miRNA]]
  # List of plots
  plots <- list()
  # Iterate stresses of a microRNA
  for (stress in names(miRNAs_list)){
    miRNAs_stress <- miRNAs_list[[stress]]
    
    # Create a dataframe for save a family of microRNA sequences expression profile through time
    df_expression_plot <- data.frame()
	  
    # Iterate sequences of a microRNA family
    for(i in 1:length(miRNAs_stress)){
      sequence <- miRNAs_stress[[i]]
      df_stress <- dataframes_micro_raw[[stress]]
      row <- data.frame(seq = sequence, row.names = paste0("seq",i), stringsAsFactors = FALSE)
      
      # Obtain the expression profile of the sequence
      for(time in 1:length(df_stress)){
        df_mirna_raw <- df_stress[[time]]
        name_lfc <- paste0("LFC_T",time)
        name_padj <- paste0("Padj_T",time)
        row[[name_lfc]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == sequence, "Shrunkenlog2FoldChange"])
        row[[name_padj]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == sequence, "padj"])
      }  
      
      # Change NA in padj to 1
      for (col in names(row)) {
        if (startsWith(col, "Padj")) {
          row[[col]][is.na(row[[col]])] <- 1
        }
      }
	    
      # If sequence has not information in some time, it is invalid
      if(!any(is.na(row))){
        
        # Add row to the final dataframe 
        df_expression_plot <- rbind(df_expression_plot, row)
        
        # Save valid microRNA 
        row_valid <- data.frame(ID = paste0("seq",i), seq=sequence, microRNA = miRNA, stress = stress, Type = "Valid",stringsAsFactors = FALSE )
        micro_summary_plot <- rbind(micro_summary_plot, row_valid) 
      }
      
      else{
        # Save invalid microRNA 
        row_invalid <- data.frame(ID = paste0("seq",i), seq=sequence, microRNA = miRNA,stress = stress, Type = "Invalid", stringsAsFactors = FALSE )
        micro_summary_plot <- rbind(micro_summary_plot, row_invalid) 
      }
    	# Add all the sequences to the final table without filter
    	row$miRNA <- miRNA
    	row$Stress <- stress
    	micro_summary_table <- rbind(micro_summary_table, row)    
    }

    if(length(df_expression_plot) != 0){
      # Change dataframe to long dataframe to do the graph
      df_long <- df_expression_plot %>%
        pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
        pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
        filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
        dplyr::select(-time_padj)

      # Create a new column for significance
      df_long$significance <- ifelse(df_long$Padj < 0.05, "Significative", "No Significative")
      
      # Sort 'seq' column
      df_long$seq <- factor(df_long$seq, levels = sort(unique(df_long$seq)))
      
      # Calculate the maximun absolute value of LFC to use it as limits in aes y
      max_abs_lfc <- max(abs(df_long$LFC))
      y_limits <- c(-max_abs_lfc - 0.5, max_abs_lfc + 0.5)
      
      # Create plot
      p <- ggplot(df_long, aes(x = time, y = LFC, group = seq, color = seq)) +
        geom_line() +
        geom_point(aes(shape = significance),size=3) +
        labs(x = "Time", y = "Log Fold Change (LFC)",color = miRNA) +
        theme_bw() +
        theme(strip.text = element_text(size = 15),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 20),
              plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15)) +
        scale_x_discrete(labels = function(x) str_replace(x, "LFC_T", "T")) +
        scale_shape_manual(name = "Statistical significance",
                           values = c("Significative" = 16, "No Significative" = 1),
                           labels = c("Significative" = "Significative", "No Significative" = "No Significative"))+
        scale_y_continuous(limits = y_limits) +
        guides(color = guide_legend(order = 1, title = miRNA)) +
        ggtitle(toupper(stress))

        # Divide the plot space into legend and graph
        p_no_legend <- p + theme(legend.position = "none")
        legend <- get_legend(p)
        combined_plot <- plot_grid(p_no_legend, legend, ncol = 2, rel_widths = c(2, 1))
      
        plots[[stress]] <- combined_plot 
    }
  }
  if(length(plots) > 0){
    final_plot <- grid.arrange(grobs = plots, ncol = 1)
    # Create and save plot with all stresses
    ggsave(plot = final_plot,filename = paste0(path_out_analysis1,"/Expression_profile_",miRNA,".svg"),
         width = 30, height = 22, dpi = 300,bg = "white")
  }
  else{
    print(paste(miRNA, "has not valid sequence in any stress"))
  }
}
# Select only the significative sequences
summary_long <- micro_summary_table %>%
  pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
  pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
  filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
  dplyr::select(-time_padj)

summary_long_filt <- summary_long[summary_long$Padj < 0.05,]

# Calculate % of significative sequences which follow the same patron for each time
expression_summary <- summary_long_filt %>%
  group_by(time,Stress,miRNA) %>%
  summarise(neg_count = sum(LFC < 0), 
            pos_count = sum(LFC > 0))	

# Save sequence table
write.table(expression_summary,paste0(path_out_logs,"/Summary_table_microRNAs.tsv"), sep='\t', col.names = TRUE, row.names = FALSE,quote=FALSE)

################################ Analysis 2 ####################################

# Extract microRNA sequences for stress
miRNAs_annot <- list()
for(stress in names(dataframes_micro_sig)){
  micro_sig_stress <- dataframes_micro_sig[[stress]]
  micro_sig_stress_table <- do.call(rbind, micro_sig_stress)
  miRNAs_annot[[stress]] <- unique(micro_sig_stress_table$general_annot)
}

# Select common microRNA
miRNA_freq <- table(unlist(miRNAs_annot))

# Filter the microRNAs that appear DE in at least 3 of the 4 stresses
miRNAs_common <- names(miRNA_freq[miRNA_freq >= 3])

# List of sequences with high BaseMean of each general miRNA in each time
general_miRNAs <- list()
for(mirna in miRNAs_common){
  temporal <- data.frame()
  
  for(stress in names(dataframes_micro_sig)){
    micro_sig_stress <- dataframes_micro_sig[[stress]]
    
    for(time in 1:length(micro_sig_stress)){
      df_time <- micro_sig_stress[[time]]
      df_time_mirna <- df_time[df_time$general_annot == mirna,]
      
      if(nrow(df_time_mirna) != 0 ){
        # Select the sequence with high baseMean at each time
        seq <- df_time_mirna[df_time_mirna$baseMean == max(df_time_mirna$baseMean),]
        
        if( stress == "drought"){
        time <- time + 1
        }
        row <- data.frame(baseMean= seq$baseMean, sequence = seq$seq, stress = stress, time = time)
        temporal <- rbind(row, temporal)
      }
    }
  }
  df_filtered <- temporal %>%
    group_by(time) %>%
    filter(baseMean == max(baseMean)) %>%
    ungroup()
  general_miRNAs[[mirna]] <- unique(df_filtered$sequence)
}

# Extract microRNA targets for stress
lista_targets <- list()
for(mirna in miRNAs_common){
  targets_tab <- targets[targets$microRNA == mirna,]
  targets_list <- targets_tab$gene
  lista_targets[[mirna]] <- targets_list
}

# Create a correlation table
correlation_table <- data.frame()
# Graph each sequence with each one of their targets
for(mirna in names(general_miRNAs)){
  mirna_list <- general_miRNAs[[mirna]]
  
  for(i in 1:length(mirna_list)){
    secuencia <- mirna_list[[i]]
    df_expression <- data.frame()
    
    for(stress in names(dataframes_micro_raw)){
      df_stress <- dataframes_micro_raw[[stress]]
      row <- data.frame(seq = secuencia, row.names = stress, stringsAsFactors = FALSE)
      
      for(time in 1:length(df_stress)){
        df_mirna_raw <- df_stress[[time]]
        name_lfc <- paste0("LFC_T",time)
        name_padj <- paste0("Padj_T",time)
        row[[name_lfc]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == secuencia, "Shrunkenlog2FoldChange"])
        row[[name_padj]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == secuencia, "padj"])
      }  
      
      # Add row to dataframe
      df_expression <- rbind(df_expression, row)
      
    }
    
    # Continue only if the sequence has not NA's
    if(!any(apply(df_expression,1,is.na))){
      # Convertir los rownames a una columna
      df_expression <- df_expression %>%
        rownames_to_column(var = "stress")
      
      # Tranform dataframe to long dataframe
      df_long <- df_expression %>%
        pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
        pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
        filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
        dplyr::select(-time_padj)
      
      # Create the tables to the targets
      if (length(lista_targets[[mirna]]) > 0 ){
        for(gene in unique(lista_targets[[mirna]])){
          df_expression_gene <- data.frame()
          
          for(stress in names(dataframes_trans)){
            df_stress <- dataframes_trans[[stress]]
            row <- data.frame(seq = gene, row.names = stress, stringsAsFactors = FALSE)
            
            for(time in 1:length(df_stress)){
              df_trans <- df_stress[[time]]
              name_lfc <- paste0("LFC_T",time)
              name_padj <- paste0("Padj_T",time)
              row[[name_lfc]] <- as.numeric(df_trans[df_trans$seq == gene, "Shrunkenlog2FoldChange"])
              row[[name_padj]] <- as.numeric(df_trans[df_trans$seq == gene, "padj"])
            }  
            
            # Add row
            df_expression_gene <- rbind(df_expression_gene, row)
          }
          # Only the genes which have not NA's are valids
          if(!any(apply(df_expression_gene,1,is.na))){
            
            # Rownames to column
            df_expression_gene <- df_expression_gene %>%
              rownames_to_column(var = "stress")
            
            # Dataframe to long dataframe
            df_long_gene <- df_expression_gene %>%
              pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
              pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
              filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
              dplyr::select(-time_padj)
            
            df_long_complete <- rbind(df_long, df_long_gene)
            
            # Create a new column to describe significance
            df_long_complete$significance <- ifelse(df_long_complete$Padj < 0.05, "Significativo", "No Significativo")
            
            # Create a new column which show if is a gene or microRNA row
            df_long_complete <- df_long_complete %>%
              mutate(shape_group = ifelse(str_detect(seq, "MELO"), "Gene", "microRNA"))
            
            # Create plot
            # Define colours
            color_values <- ifelse(unique(df_long_complete$seq) %in% df_long_complete$seq[str_detect(df_long_complete$seq, "MELO")], "cadetblue2", "#FF6347")
            names(color_values) <- unique(df_long_complete$seq)

            # Calculate the aes limits
            max_abs_lfc <- max(abs(df_long_complete$LFC))
            y_limits <- c(-max_abs_lfc - 0.5, max_abs_lfc + 0.5)
            
            # Plot
            p <- ggplot(df_long_complete, aes(x = time, y = LFC, group = seq, color = seq, shape = interaction(shape_group, significance), fill = interaction(shape_group, significance))) +
              geom_line() +
              geom_point(size = 3) +
              facet_wrap(~ stress, scales = "free_y",
                         labeller = labeller(stress = c("cold" = "COLD", "drought" = "DROUGHT", "mon" = "MONOSPORASCUS", "sd" = "SHORTDAY")), nrow = 1) +
              labs(x = "Time", y = "Log Fold Change (LFC)", color = mirna, fill = "Sequence type", fill = "Significance") +
              theme_bw() +
              theme(strip.text = element_text(size = 10)) +
              scale_x_discrete(labels = function(x) str_replace(x, "LFC_T", "T")) +
              scale_shape_manual(name = "Sequence & Significance",
                                 values = c("microRNA.Significativo" = 24, "microRNA.No Significativo" = 2, 
                                            "Gene.Significativo" = 21, "Gene.No Significativo" = 1), 
                                 labels = c("microRNA.Significativo" = "Significative microRNA", "microRNA.No Significativo" = "No significative microRNA",
                                            "Gene.Significativo" = "Significative Gene", "Gene.No Significativo" = "No significative Gene")) +
              scale_fill_manual(name = "Sequence & Significance",
                                values = c("microRNA.Significativo" = "#FF6347", "microRNA.No Significativo" = "white",
                                           "Gene.Significativo" = "cadetblue3", "Gene.No Significativo" = "white"),
                                labels = c("microRNA.Significativo" = "Significative microRNA", "microRNA.No Significativo" = "No significative microRNA",
                                           "Gene.Significativo" = "Significative Gene", "Gene.No Significativo" = "No significative Gene")) +
              scale_color_manual(name = "miRNA",
                                 values = color_values) +
              scale_y_continuous(limits = y_limits) +
              guides(color = guide_legend(order = 1, title = mirna))
            
            ggsave(plot = p,filename = paste0(path_out_analysis2,"/Expression_profile_",mirna,"_",i,"_",gene,".svg"),
                   width = 20, height = 7, dpi = 300,bg = "white")
            
            # Create row of correlation matrix only with lfc
	          df_long_complete_sig <- df_long_complete[df_long_complete$significance == "Significativo",]
            lfc_row_gene <- df_long_complete_sig[df_long_complete_sig$shape_group == "Gene", c("LFC","time", "stress")]
            colnames(lfc_row_gene) <- c("LFC_gene", "time","stress")
            lfc_row_micro <- df_long_complete_sig[df_long_complete_sig$shape_group == "microRNA", c("LFC","time", "stress")]
            colnames(lfc_row_micro) <- c("LFC_micro", "time","stress")
            lfc_row <- merge(lfc_row_gene, lfc_row_micro, by = c("time","stress"))
            
            if (nrow(lfc_row) != 0){
              lfc_row$microRNA <- mirna
              lfc_row$Gene <- gene
              # Agregar la fila a la tabla de correlación
              correlation_table <- rbind(correlation_table, lfc_row)     
            }
          } 
          
          else{
          print(paste0("The gene ", gene, " have NA's"))
            
          }
        }
      }
      
      else {
        print(paste0("The sequence ", secuencia, " of ", mirna, " have not valid targets"))
        
      }
    }
    
    else{
      print(paste0("The sequence ", secuencia, ", number ", i, " of ", mirna, " have NA's"))
      
    }
  }
}

print("Claculating the correlation...")

write.table(correlation_table,paste0(path_out_analysis2,"/Correlation_table.tsv"), row.names= FALSE, col.names = TRUE, quote= FALSE)

# Calculate the correlation to a No normal data distribution
correlation_table$LFC_gene_jitter <- jitter(correlation_table$LFC_gene)
correlation_table$LFC_micro_jitter <- jitter(correlation_table$LFC_micro)
cor_spearman <- cor.test(correlation_table$LFC_gene_jitter, correlation_table$LFC_micro_jitter, method = "spearman")

# Clculate the max malue to the axes
max_abs_x <- max(abs(correlation_table$LFC_micro), na.rm = TRUE) +0.5
max_abs_y <- max(abs(correlation_table$LFC_gene), na.rm = TRUE) + 0.5

# Extract the Spearman correlation coeficient
spearman_coefficient <- cor_spearman$estimate
spearman_p <- cor_spearman$p.value

correlation_plot <- ggplot(correlation_table, aes(x = LFC_micro, y = LFC_gene, color = stress)) +
  geom_point(size = 2) +
  xlim(-max_abs_x, max_abs_x) +  # Ajusta los límites del eje x
  ylim(-max_abs_y, max_abs_y) +  # Ajusta los límites del eje y
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Línea vertical en x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea horizontal en y = 0
  annotate("text", 
           label = paste("Spearman's correlation:", round(spearman_coefficient, digits = 4), 
                         "\np-value:", round(spearman_p, digits = 10)), 
           x = 7, y = 6, size = 4, color = "black") +  # Añadir un comentario
  theme_classic() +  # Estilo de tema minimalista
  labs(x = "microRNA LFC", y = "Genes LFC", color = "Stress")  +  # Etiquetas de los ejes y la leyenda
  theme(
    plot.title = element_text(size = 16),      # Tamaño y estilo del título del gráfico
    axis.title.x = element_text(size = 14, face = "bold"),                    # Tamaño del título del eje X
    axis.title.y = element_text(size = 14, face = "bold"),                    # Tamaño del título del eje Y
    axis.text.x = element_text(size = 12),                     # Tamaño de las etiquetas del eje X
    axis.text.y = element_text(size = 12),                     # Tamaño de las etiquetas del eje Y
    legend.title = element_text(size = 14),     # Tamaño y estilo del título de la leyenda
    legend.text = element_text(size = 12),                     # Tamaño del texto de la leyenda
    legend.position = "top"                                    # Colocar la leyenda en la parte superior
  )

ggsave(paste0(path_out_network,"/CORRELATION_PLOT.png"), 
       plot = correlation_plot, width = 10, height = 10, bg = "white")

########################## HEATMAP AND NETWORK #################################
print("Generating Network...")
# Extract the time information from the column
correlation_table$time <- gsub(".*_", "", correlation_table$time)

# Filter rows where LFC_gene and LFC_micro have opposite signs
correlation_table_filt <- correlation_table %>%
  filter((LFC_gene > 0 & LFC_micro < 0) | (LFC_gene < 0 & LFC_micro > 0))

# Create the edges table with the count of each unique miRNA-Gene combination
edge_table <- correlation_table_filt %>%
  group_by(microRNA, Gene) %>%
  summarise(count = n(), .groups = 'drop')

# Rename columns to match edge list format
colnames(edge_table) <- c("target", "source", "weight")

# Save the edges table to a file
write.table(edge_table, 
            paste0(path_out_network,"/Edge_table.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Create the nodes table for genes
genes <- unique(correlation_table_filt$Gene)
table_genes <- data.frame()
for (gene in genes){
  n <- 0
  for (condition in names(dataframes_trans)){
    tables <- dataframes_trans[[condition]]
    for (i in seq_along(tables)){
      table <- tables[[i]]
      if (gene %in% table$seq){
        if (table$padj[table$seq == gene] < 0.05){
          n <- n + 1
        }
      }
    }
  }
  row <- data.frame(id= gene, score = n, group = "gene")
  table_genes <- rbind(table_genes,row)
}

# Create the nodes table for microRNAs
microRNA <- unique(correlation_table_filt$microRNA)
table_micro <- data.frame()
for (micro in microRNA){
  n <- 0
  for (condition in names(dataframes_micro_sig)){
    tables <- dataframes_micro_sig[[condition]]
    for (i in seq_along(tables)){
      table <- tables[[i]]
      if (micro %in% table$general_annot){
        if (any(table$padj[table$general_annot == micro] < 0.05)){
          n <- n + 1
        }
      }
    }
  }
  row <- data.frame(id= micro, score = n, group = "microRNA")
  table_micro <- rbind(table_micro,row)
}

# Combine gene and microRNA node tables
nodes_table <- rbind(table_genes, table_micro)

# Save the nodes table to a file
write.table(nodes_table, 
            paste0(path_out_network,"/Nodes_table.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
            
# Group Go
genes <- annotation_table$gene
ggo <- groupGO(gene = table_genes$id,
                         OrgDb =org.CMelo.eg.db,
                         keyType = "GID",
                         ont      = "MF",
                         level    = 3)
results <- ggo@result
results <- results[results$Count != 0, ]
# Save the Go terms table to a file
write.table(results, 
            paste0(path_out_network,"/GO_terms_MF_3.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ggo <- groupGO(gene = table_genes$id,
                         OrgDb =org.CMelo.eg.db,
                         keyType = "GID",
                         ont      = "BP",
                         level    = 3)
results <- ggo@result
results <- results[results$Count != 0, ]
# Save the Go terms table to a file
write.table(results, 
            paste0(path_out_network,"/GO_terms_BP_3.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

print("Generating heatmap...")
############## microRNA Processing ##############
# Select relevant columns for microRNA
micro_table <- correlation_table_filt[, c("time", "stress", "LFC_micro", "microRNA")]

# Remove duplicate rows
micro_table_uniq <- micro_table %>% distinct()

# Combine "stress" and "time" into a single column "stress_time"
micro_table_uniq <- micro_table_uniq %>%
  unite("stress_time", stress, time, sep = "_")

# Select the row with the highest absolute LFC_micro per microRNA and stress_time
micro_table_uniq <- micro_table_uniq %>%
  group_by(microRNA, stress_time) %>%
  dplyr::slice(which.max(abs(LFC_micro))) %>%
  ungroup()

# Define stress conditions and time points
stresses <- unique(correlation_table_filt$stress)
times <-  unique(correlation_table_filt$time)
mirnas <- unique(micro_table_uniq$microRNA)

# Create an empty dataframe to store microRNA matrix
matrix <- data.frame()

# Iterate over each microRNA to populate the matrix
for (mirna in mirnas) {
  lfc_table <- micro_table_uniq[micro_table_uniq$microRNA == mirna, ]
  row <- data.frame(matrix(0, nrow = 1, ncol = length(stresses) * length(times)))
  colnames(row) <- c(outer(stresses, times, paste, sep = "_"))
  
  for (time in times) {
    for (stress in stresses) {
      name <- paste0(stress, "_", time)
      row[[name]] <- ifelse(name %in% lfc_table$stress_time, as.numeric(lfc_table[lfc_table$stress_time == name, "LFC_micro"]), NA)
    }
  }
  rownames(row) <- mirna
  matrix <- rbind(matrix, row)
}

############## Gene Processing ##############
# Select relevant columns for genes
gene_table <- correlation_table_filt[, c("time", "stress", "LFC_gene", "Gene")]

gene_table_uniq <- gene_table %>% distinct()

gene_table_uniq <- gene_table_uniq %>%
  unite("stress_time", stress, time, sep = "_")

gene_table_uniq <- gene_table_uniq %>%
  group_by(Gene, stress_time) %>%
  dplyr::slice(which.max(abs(LFC_gene))) %>%
  ungroup()

# Define unique genes
genes <- unique(gene_table_uniq$Gene)

# Create an empty dataframe to store gene matrix
matrix_2 <- data.frame()

# Iterate over each gene to populate the matrix
for (gene in genes) {
  lfc_table <- gene_table_uniq[gene_table_uniq$Gene == gene, ]
  row <- data.frame(matrix(0, nrow = 1, ncol = length(stresses) * length(times)))
  colnames(row) <- c(outer(stresses, times, paste, sep = "_"))
  print(head(lfc_table))
  for (time in times) {
    for (stress in stresses) {
      name <- paste0(stress, "_", time)
      row[[name]] <- ifelse(name %in% lfc_table$stress_time, as.numeric(lfc_table[lfc_table$stress_time == name, "LFC_gene"]), NA)
      print(head(row))
    }
  }
  rownames(row) <- gene
  matrix_2 <- rbind(matrix_2, row)
}

############## Annotation Tables ##############
# Create annotation table for microRNA
annotation_microRNA <- data.frame(row.names = rownames(matrix), microRNA = rownames(matrix))
annotation_microRNA <- annotation_microRNA %>%
  mutate(microRNA = ifelse(microRNA %in% c("miR156", "miR157"), "miR156-miR157", microRNA))

# Load target annotation
annotation_microRNA[["target"]] <- rownames(annotation_microRNA)
annotation_gene <- merge(annotation_microRNA, edge_table, by = "target")
annotation_gene <- annotation_gene[, c("microRNA", "source")]
annotation_gene <- annotation_gene %>% distinct()
rownames(annotation_gene) <- annotation_gene$source
colnames(annotation_gene) <- c("microRNA", "other")
colnames(annotation_microRNA) <- c("microRNA", "other")
annotation_table <- rbind(annotation_gene, annotation_microRNA)

############## Heatmap Construction ##############
# Combine microRNA and gene matrices
HM_matrix <- rbind(matrix, matrix_2)
HM_matrix <- as.matrix(HM_matrix)

# Define row order
rownames_order <- c("miR164", "miR319", "miR166", "miR157", "miR156",
                    "miR396", "miR398", "miR408", "MELO3C017185", "MELO3C007121", "MELO3C002754",
                    "MELO3C016092", "MELO3C007078", "MELO3C017245",
                    "MELO3C007656", "MELO3C009159", "MELO3C022680", "MELO3C016781",
                    "MELO3C015374", "MELO3C008424", "MELO3C027302")

# Order columns and rows
HM_matrix <- HM_matrix[, order(colnames(HM_matrix))]
HM_matrix <- HM_matrix[match(rownames_order, rownames(HM_matrix)), ]
annotation_table <- annotation_table[match(rownames_order, rownames(annotation_table)), ]

# Define color function
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun(seq(-1, 1))

# Define microRNA annotation colors
microRNA_annot_colors <- c("miR156-miR157" = "#FF7F00", "miR319" = "#FFFF32",
                           "miR166" = "#32FF00", "miR396" = "#A5EDFF", "miR164" = "#CCBFFF",
                           "miR408" = "#654CFF", "miR398" = "#E51932")

# Create row annotation
annotation_row <- rowAnnotation(microRNA = annotation_table$microRNA, col = list(microRNA = microRNA_annot_colors))

row_group <- ifelse(grepl("^miR", rownames_order), "miR", "MELO")

HM_matrix_filtrada <- HM_matrix[, colSums(!is.na(HM_matrix)) > 0]
print(head(HM_matrix_filtrada))
heat <- Heatmap(HM_matrix_filtrada, rect_gp = gpar(col = "white", lwd = 0.5),
                col = col_fun, right_annotation = annotation_row, 
                column_names_rot = 45, column_names_gp = gpar(fontsize = 16),
                row_names_gp = gpar(fontsize = 16),
                column_labels = c("Cold T1","Cold T2","Cold T3","Drought T3",
                                  "Monosporascus T2","ShortDay T2","ShortDay T3"),
                row_split = row_group,
                gap = unit(2, "mm"),
                row_title = c("Genes","microRNA" ),
                cluster_rows = FALSE, cluster_columns = FALSE)

svg(paste0(path_out_network,"/Heatmap.svg"), width = 17, height = 10)
heat
dev.off()

############################### Analysis 3 #####################################

miRNA_diff <- list(miR156=list(),miR159=list(),miR166=list(),miR319=list(),miR396=list())

for(stress in names(dataframes_micro_sig)){
  micro_sig_stress <- dataframes_micro_sig[[stress]]
  for(time in 1:length(micro_sig_stress)){
    df_time <- micro_sig_stress[[time]]
    for(mirna in names(miRNA_diff)){
      lista_seq <- df_time[df_time$general_annot == mirna,]
      lista_seq <- lista_seq$seq
      miRNA_diff[[mirna]] <- as.character(unique(c(miRNA_diff[[mirna]], lista_seq)))
    }
  }
}

# Graph each sequence with each one of their targets
for(mirna in names(miRNA_diff)){
  mirna_list <- miRNA_diff[[mirna]]
  
  for(i in 1:length(mirna_list)){
    secuencia <- mirna_list[[i]]
    df_expression <- data.frame()
    
    for(stress in names(dataframes_micro_raw)){
      df_stress <- dataframes_micro_raw[[stress]]
      row <- data.frame(seq = secuencia, row.names = stress, stringsAsFactors = FALSE)
      
      for(time in 1:length(df_stress)){
        df_mirna_raw <- df_stress[[time]]
        name_lfc <- paste0("LFC_T",time)
        name_padj <- paste0("Padj_T",time)
        row[[name_lfc]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == secuencia, "Shrunkenlog2FoldChange"])
        row[[name_padj]] <- as.numeric(df_mirna_raw[df_mirna_raw$seq == secuencia, "padj"])
      }  
      
      # Add row to dataframe
      df_expression <- rbind(df_expression, row)
      
    }
    
    # Continue only if the sequence has not NA's
    if(!any(apply(df_expression,1,is.na))){
      # Convertir los rownames a una columna
      df_expression <- df_expression %>%
        rownames_to_column(var = "stress")
      
      # Tranform dataframe to long dataframe
      df_long <- df_expression %>%
        pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
        pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
        filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
        dplyr::select(-time_padj)
      
      # Create the tables to the targets
      if (length(lista_targets[[mirna]]) > 0 ){
        for(gene in unique(lista_targets[[mirna]])){
          df_expression_gene <- data.frame()
          
          for(stress in names(dataframes_trans)){
            df_stress <- dataframes_trans[[stress]]
            row <- data.frame(seq = gene, row.names = stress, stringsAsFactors = FALSE)
            
            for(time in 1:length(df_stress)){
              df_trans <- df_stress[[time]]
              name_lfc <- paste0("LFC_T",time)
              name_padj <- paste0("Padj_T",time)
              row[[name_lfc]] <- as.numeric(df_trans[df_trans$seq == gene, "Shrunkenlog2FoldChange"])
              row[[name_padj]] <- as.numeric(df_trans[df_trans$seq == gene, "padj"])
            }  
            
            # Add row
            df_expression_gene <- rbind(df_expression_gene, row)
            
          }
          # Only the genes which have not NA's are valids
          if(!any(apply(df_expression_gene,1,is.na))){
            
            # Rownames to column
            df_expression_gene <- df_expression_gene %>%
              rownames_to_column(var = "stress")
            
            # Dataframe to long dataframe
            df_long_gene <- df_expression_gene %>%
              pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
              pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
              filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
              dplyr::select(-time_padj)
            
            df_long_complete <- rbind(df_long, df_long_gene)
            
            # Create a new column to describe significance
            df_long_complete$significance <- ifelse(df_long_complete$Padj < 0.05, "Significativo", "No Significativo")
            
            # Create a new column which show if is a gene or microRNA row
            df_long_complete <- df_long_complete %>%
              mutate(shape_group = ifelse(str_detect(seq, "MELO"), "Gene", "microRNA"))
            
            # Create plot
            # Define colours
            color_values <- ifelse(unique(df_long_complete$seq) %in% df_long_complete$seq[str_detect(df_long_complete$seq, "MELO")], "cadetblue2", "#FF6347")
            names(color_values) <- unique(df_long_complete$seq)
            
            # Calculate the aes limits
            max_abs_lfc <- max(abs(df_long_complete$LFC))
            y_limits <- c(-max_abs_lfc - 0.5, max_abs_lfc + 0.5)
            
            # Plot
            p <- ggplot(df_long_complete, aes(x = time, y = LFC, group = seq, color = seq, shape = interaction(shape_group, significance), fill = interaction(shape_group, significance))) +
              geom_line() +
              geom_point(size = 3) +
              facet_wrap(~ stress, scales = "free_y",
                         labeller = labeller(stress = c("cold" = "COLD", "drought" = "DROUGHT", "mon" = "MONOSPORASCUS", "sd" = "SHORTDAY")), nrow = 1) +
              labs(x = "Time", y = "Log Fold Change (LFC)", color = mirna, fill = "Sequence type", fill = "Significance") +
              theme_bw() +
              theme(strip.text = element_text(size = 10)) +
              scale_x_discrete(labels = function(x) str_replace(x, "LFC_T", "T")) +
              scale_shape_manual(name = "Sequence & Significance",
                                 values = c("microRNA.Significativo" = 24, "microRNA.No Significativo" = 2, 
                                            "Gene.Significativo" = 21, "Gene.No Significativo" = 1), 
                                 labels = c("microRNA.Significativo" = "Significative microRNA", "microRNA.No Significativo" = "No significative microRNA",
                                            "Gene.Significativo" = "Significative Gene", "Gene.No Significativo" = "No significative Gene")) +
              scale_fill_manual(name = "Sequence & Significance",
                                values = c("microRNA.Significativo" = "#FF6347", "microRNA.No Significativo" = "white",
                                           "Gene.Significativo" = "cadetblue3", "Gene.No Significativo" = "white"),
                                labels = c("microRNA.Significativo" = "Significative microRNA", "microRNA.No Significativo" = "No significative microRNA",
                                           "Gene.Significativo" = "Significative Gene", "Gene.No Significativo" = "No significative Gene")) +
              scale_color_manual(name = "miRNA",
                                 values = color_values) +
              scale_y_continuous(limits = y_limits) +
              guides(color = guide_legend(order = 1, title = mirna))
            
            ggsave(plot = p,filename = paste0(path_out_analysis3,"/Expression_profile_",mirna,"_",i,"_",gene,".svg"),
                   width = 20, height = 7, dpi = 300,bg = "white")
          } 
          
          else{
            print(paste0("The gene ", gene, " have NA's"))
            
          }
        }
      }
      
      else {
        print(paste0("The sequence ", secuencia, " of ", mirna, " have not valid targets"))
        
      }
    }
    
    else{
      print(paste0("The sequence ", secuencia, ", number ", i, " of ", mirna, " have NA's"))
      
    }
  }
}

