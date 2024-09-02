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
# and is graphicated with each one of their targets profile.
# - Analysis 3: Select only the microRNAs that have a different profile that
# their family and and graphicated with each one of their targets profile.

rm(list = ls())

################################# Libraries ####################################
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))
suppressMessages(library("argparse"))
suppressMessages(library(svglite))
suppressMessages(library("cowplot"))
suppressMessages(library("plotly"))
suppressMessages(library(htmlwidgets))

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

# Create directories if they do not exist
dir.create(path_out_analysis1 , recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_analysis2, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_analysis3, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_logs, recursive = TRUE, showWarnings = FALSE)

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

############################ Analysis 1 ########################################

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
micro_summary <- data.frame()

# Iterate microRNAs
for(miRNA in names(miRNAs)){
  miRNAs_list <- miRNAs[[miRNA]]
  # List of plots
  plots <- list()
  # Iterate stresses of a microRNA
  for (stress in names(miRNAs_list)){
  miRNAs_stress <- miRNAs_list[[stress]]
    
    # Create a dataframe for save a family of microRNA sequences expression profile through time
    df_expression <- data.frame()
    
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
        df_expression <- rbind(df_expression, row)
        
        # Save valid microRNA 
        row_valid <- data.frame(ID = paste0("seq",i), seq=sequence, microRNA = miRNA, stress = stress, Type = "Valid", stringsAsFactors = FALSE )
        micro_summary <- rbind(micro_summary, row_valid) 
      }
      
      else{
        # Save invalid microRNA 
        row_invalid <- data.frame(ID = paste0("seq",i), seq=sequence, microRNA = miRNA,stress = stress, Type = "Invalid", stringsAsFactors = FALSE )
        micro_summary <- rbind(micro_summary, row_invalid) 
      }
    }
    
    if(length(df_expression) != 0){
      # Change dataframe to long dataframe to do the graph
      df_long <- df_expression %>%
        pivot_longer(cols = starts_with("LFC_T"), names_to = "time", values_to = "LFC") %>%
        pivot_longer(cols = starts_with("Padj_T"), names_to = "time_padj", values_to = "Padj") %>%
        filter(str_sub(time, 5) == str_sub(time_padj, 6)) %>%
        select(-time_padj)
      
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

# Save sequence table
  write.table(micro_summary,paste0(path_out_logs,"/Summary_table_microRNAs.tsv"), sep='\t', col.names = TRUE, row.names = FALSE)

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
        # Select the microRNA with high baseMean at each time
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
        select(-time_padj)
      
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
              select(-time_padj)
            
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
            lfc_row_gene <- df_long_complete[df_long_complete$shape_group == "Gene", c("LFC","time", "stress")]
            colnames(lfc_row_gene) <- c("LFC_gene", "time","stress")
            lfc_row_micro <- df_long_complete[df_long_complete$shape_group == "microRNA", c("LFC","time", "stress")]
            colnames(lfc_row_micro) <- c("LFC_micro", "time","stress")
            lfc_row <- merge(lfc_row_gene, lfc_row_micro, by = c("time","stress"))
            lfc_row$microRNA <- mirna
            lfc_row$Gene <- gene
            
            # Save the row only if bpth LFC are higher than 0.5
            for (i in 1:nrow(lfc_row)) {
              if ((abs(lfc_row$LFC_gene[i]) >= 0.5) & (abs(lfc_row$LFC_micro[i]) >= 0.5)) {
                # Agregar la fila a la tabla de correlación
                correlation_table <- rbind(correlation_table, lfc_row[i, ])
              }
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

write.table(correlation_table,paste0(path_out_analysis2,"/Correlation_table.txt"))
# Calculate the correlation to a No normal data distribution
cor_spearman <- cor.test(correlation_table$LFC_gene, correlation_table$LFC_micro, method = "spearman")
print(cor_spearman)

cor_kendall <- cor.test(correlation_table$LFC_gene, correlation_table$LFC_micro, method = "kendall")
print(cor_kendall)

# Create the plot
p <- plot_ly(correlation_table, x = ~LFC_gene, y = ~LFC_micro,
               text = ~paste("microRNA: ", microRNA, '<br>Gene:', Gene),
               color = ~stress,
               type = 'scatter',
               mode = 'markers')

# Extract the Spearman correlation coeficient
spearman_coefficient <- cor_spearman$estimate

# Add anotation to the plot
p <- p %>% layout(
  annotations = list(
    x = 4,  # Coordenada x para la anotación
    y = 10,  # Coordenada y para la anotación
    text = paste("Correlación Spearman:", round(spearman_coefficient, 2)),  # Texto de la anotación
    showarrow = FALSE,  # Ocultar la flecha
    xref = "x",  # Referencia de la coordenada x
    yref = "y",  # Referencia de la coordenada y
    xanchor = 'left',  # Alineación horizontal del texto
    yanchor = 'bottom'  # Alineación vertical del texto
  )
)

# Save the plot as HTML file
htmlwidgets::saveWidget(p, paste0(path_out_analysis2,"/Correlation_dotplot.html"), selfcontained = TRUE)

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
        select(-time_padj)
      
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
              select(-time_padj)
            
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

