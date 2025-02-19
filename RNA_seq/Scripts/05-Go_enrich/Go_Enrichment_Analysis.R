################################################################################
##                                                                            ##
##  GO_Enrichment_Analysis.R                                                  ##
##                                                                            ##
##  1. Differential expression data processing                                ##
##                                                                            ##
##  This script reads differential expression analysis results from various   ##
##  stress conditions and time points. It extracts significantly upregulated  ##
##  and downregulated genes for each condition.                               ##
##                                                                            ##
##  2. Gene Ontology (GO) Enrichment Analysis                                 ##
##                                                                            ##
##  Using the extracted differentially expressed genes, the script performs   ##
##  GO enrichment analysis using the clusterProfiler package. The analysis    ##
##  considers biological processes (BP) and filters results based on an       ##
##  adjusted p-value threshold.                                               ##
##                                                                            ##
##  3. Visualization and Results Export                                       ##
##                                                                            ##
##  The script generates bar plots of enriched GO terms for both upregulated  ##
##  and downregulated genes. Additionally, it creates upset plots to          ##
##  visualize shared GO terms across conditions. The results are saved as     ##
##  tables and plots in specified output directories.                         ##
##                                                                            ##
##  Author: Marta Núñez Salvador                                              ##
##  Version: 1.0                                                              ##
##  Date: 18/02/2025                                                          ##
##                                                                            ##
################################################################################

suppressMessages(library(org.CMelo.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library("enrichplot"))
suppressMessages(library(readr))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexUpset))
suppressMessages(library(tidyr))
suppressMessages(library(argparse))

################################## FUNCTIONS ###################################

#' Get the command line arguments
#' This function parses the command line arguments entered into the program.
#'
#' @return List with the argument values
get_arguments <- function() {
  parser <- ArgumentParser(prog = 'GO_Enrichment_Analysis.R',
                           description = 'This script performs a Gene Ontology (GO) enrichment analysis on differentially expressed genes, 
                           generates plots, and exports results.')
  
  required <- parser$add_argument_group('required arguments')
  
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'Path to the directory containing transcript count files.',
                        required = TRUE)
  
  required$add_argument('-a', '--annotation',
                        type = 'character',
                        help = 'Path to the annotation file.',
                        required = TRUE)
  
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Path to the output directory where results will be saved.',
                        required = TRUE)
  
  args <- parser$parse_args()
  
  if (!dir.exists(args$input)) {
    stop('Error: The transcript directory does not exist.')
  }
  
  if (!file.exists(args$annotation)) {
    stop('Error: The annotation file does not exist.')
  }
  
  return(args)
}

##################################### MAIN #####################################

# Get program arguments
args <- get_arguments()

# Save the arguments in variables
transcripts_path <- args$input
annotation_file <- args$annotation
output_path <- args$output

# Create output directory if it does not exist
output_path_table<- paste0(output_path, "/01-Tables")
output_path_plot<- paste0(output_path, "/02-Plots")
output_path_common<- paste0(output_path, "/03-Common")

dir.create(output_path_table, recursive = TRUE, showWarnings = FALSE)
dir.create(output_path_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_path_common, recursive = TRUE, showWarnings = FALSE)


# Load annotation file
Anotation_file <- read.table(annotation_file, sep = '\t', header= FALSE, quote = "")
colnames(Anotation_file) <- c('Gene','Desc.')

# Load transcript files
transcripts_files <- list.files(transcripts_path, full.names = TRUE)
sig_names <- sapply(strsplit(basename(transcripts_files), "_"), `[`, 1)
transcripts <- split(transcripts_files, sig_names)

dataframes_transcritos <- list()
for(stress in names(transcripts)){
  tablas_trans <- lapply(transcripts[[stress]],read_csv)
  dataframes_transcritos[[stress]] <- tablas_trans
}

final_table_up <- data.frame()
final_table_down <- data.frame()
#  GO Enrichment analysis of differential expressed transcripts for each stress and each time
for(estres in names(dataframes_transcritos)){
  lista_dataframes <- dataframes_transcritos[[estres]]
  for(n in 1:length(lista_dataframes)) {
      transcritos_sig <- lista_dataframes[[n]]
      
      # Separate up and down regulated transcripts
      upreg <- transcritos_sig[transcritos_sig$Shrunkenlog2FoldChange > 0,]
      upreg <- upreg$seq
      downreg <- transcritos_sig[transcritos_sig$Shrunkenlog2FoldChange < 0,]
      downreg <- downreg$seq
      
      # GO enrichment analysis
      Gene_list_universe = unique(Anotation_file$Gene)
      if(length(upreg) > 10){
        upEGO = enrichGO(gene = upreg,
                         universe = Gene_list_universe,
                         OrgDb = org.CMelo.eg.db,
                         keyType = 'GID',
                         ont = "BP",
                         minGSSize = 10,
                         maxGSSize = 500,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         #qvalueCutoff = 0.05,
                         readable = TRUE,
                         pool = FALSE)
        
        
        if(nrow(upEGO@result) != 0) {
          # Use the simplify function to reduce redundancy of enriched GO terms.
          upSimGO = simplify(upEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
                             semData = NULL)
          if (nrow(upSimGO@result) != 0){
            # Plot analysis
            # up
            png(filename=paste0(output_path_plot,'/',estres,'_T',n,"_up.png"), width = 9, height = 6, units = "in", res = 300)
            print(barplot(upSimGO, showCategory = 20) +
                    ggtitle(paste0("GO ORA of up-regulated genes")) +
                    xlab("Enriched terms") + ylab("Count")+
                    theme(axis.text.y = element_text(size = 6)))
            dev.off()
            
            # Create and down tables
            up.tab = upSimGO@result
            write.table(up.tab, file = paste0(output_path_table,'/',estres,'_T',n, "_up.txt"), sep = "\t", quote = F, 
                         row.names = F, col.names = T) 
            up_name <- paste0("up_",estres,"_T",n)
            up.tab$Estres <- paste0(estres)
            final_table_up <- rbind(final_table_up,up.tab)
          }
        }
      }
      
      if(length(downreg) > 10){
        downEGO = enrichGO(gene = downreg,
                           universe = Gene_list_universe,
                           OrgDb = org.CMelo.eg.db,
                           keyType = 'SYMBOL',
                           ont = "BP",
                           minGSSize = 10,
                           maxGSSize = 500,
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           #qvalueCutoff = 0.05,
                           readable = TRUE)
        
        if(nrow(downEGO@result) != 0) {
          # Use the simplify function to reduce redundancy of enriched GO terms.
          downSimGO = simplify(downEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang",
                               semData = NULL)
          if (nrow(downSimGO@result) != 0){
            # Plot analysis
            # down
            png(filename= paste0(output_path_plot,'/',estres,'_T',n,"_down.png"), width = 9, height = 6, units = "in", res = 300)
            print(barplot(downSimGO, showCategory = 20) +
                    ggtitle(paste0("GO ORA of down-regulated genes")) +
                    xlab("Enriched terms") + ylab("Count")+
                    theme(axis.text.y = element_text(size = 6)))
            dev.off()
            
            # Create and down tables
            dn.tab = downSimGO@result    
            write.table(dn.tab, file = paste0(output_path_table,'/',estres,'_T',n,"_down.txt"), sep = "\t", quote = F, 
                         row.names = F, col.names = T)
            down_name <- paste0("down_",estres,"_T",n)
            dn.tab$Estres <- paste0(estres)
            final_table_down <- rbind(final_table_down,dn.tab)
          }
        }
      }
  }
}



upset_up <- final_table_up %>%
  select(ID, Estres) %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Estres, values_from = Presence, values_fill = list(Presence = 0))

upset_up <- upset_up %>%
  select(ID, sort(colnames(upset_up)[-1]))

stresses = colnames(upset_up)[-1]
# transfrom in a boolean matrix
upset_up[stresses] = upset_up[stresses] == 1

stresses <- sort(stresses, decreasing = TRUE)

plot_upset_up <- upset(upset_up, stresses, name='stresses', width_ratio=0.1, height_ratio=1,sort_sets = FALSE,
                    base_annotations=list('Intersection size'=intersection_size(
                      text=list(vjust=-0.1, hjust=-0.1,angle=45, color="grey"))),
                    themes=upset_modify_themes(
                      list(
                        'main_bar' = theme(text = element_text(size = 20)),         # Tamaño de texto en la barra principal
                        'sets' = theme(text = element_text(size = 20)),             # Tamaño de texto en las etiquetas de conjuntos
                        'intersections_matrix' = theme(text = element_text(size = 15)), # Tamaño de texto en la matriz de intersecciones
                        'intersection_sizes' = theme(text = element_text(size = 20)),   # Tamaño de texto en tamaños de intersección
                        'sets_sizes' = theme(text = element_text(size = 20))             # Tamaño de texto en tamaños de conjuntos
                        
                      )
                    ))


id_down <- upset_down %>%
  rowwise() %>%
  filter(sum(c_across(-ID)) == 3) %>%
  pull(ID)

table_down <- unique(final_table_down[final_table_down$ID %in% id_down, c("ID","Description")])

id_up <- upset_up %>%
  rowwise() %>%
  filter(sum(c_across(-ID)) == 3) %>%
  pull(ID)

table_up <- unique(final_table_up[final_table_up$ID %in% id_up, c("ID","Description")])

write.table(table_up, file = paste0(output_path_common,'/common_BP_up.txt'), sep = "\t", quote = F, 
            row.names = F, col.names = T) 
write.table(table_down, file = paste0(output_path_common,'/common_BP_down.txt'), sep = "\t", quote = F, 
            row.names = F, col.names = T) 


ggsave(output_path_common,'/upset_down_plot.png', plot = plot_upset_down, width = 15, height = 10)
ggsave(output_path_common,'/upset_up_plot.png', plot = plot_upset_up, width = 15, height = 10)

####################### UPSET PLOT FOR UP-REGULATED GENES ######################
# Prepare data for the upset plot
upset_up <- final_table_up %>%
  select(ID, Estres) %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Estres, values_from = Presence, values_fill = list(Presence = 0))

# Sort columns to maintain consistent ordering
upset_up <- upset_up %>%
  select(ID, sort(colnames(upset_up)[-1]))

# Extract stress condition names
stresses = colnames(upset_up)[-1]

# Convert to boolean matrix for upset plot compatibility
upset_up[stresses] = upset_up[stresses] == 1

# Sort stress conditions in descending order for better visualization
stresses <- sort(stresses, decreasing = TRUE)

# Generate the upset plot
plot_upset_up <- upset(
  upset_up, stresses, name='stresses', width_ratio=0.1, height_ratio=1, sort_sets = FALSE,
  base_annotations=list('Intersection size' = intersection_size(
    text = list(vjust = -0.1, hjust = -0.1, angle = 45, color = "grey"))),
  themes = upset_modify_themes(
    list(
      'main_bar' = theme(text = element_text(size = 20)),   
      'sets' = theme(text = element_text(size = 20)),      
      'intersections_matrix' = theme(text = element_text(size = 15)),
      'intersection_sizes' = theme(text = element_text(size = 20)),  
      'sets_sizes' = theme(text = element_text(size = 20))          
    )
  )
)

ggsave(paste0(output_path_common, '/upset_up_plot.png'), plot = plot_upset_up, width = 15, height = 10)

# Identify common up-regulated biological processes in three conditions
id_up <- upset_up %>%
  rowwise() %>%
  filter(sum(c_across(-ID)) == 3) %>%
  pull(ID)

# Extract the unique gene IDs and their descriptions for up-regulated processes
table_up <- unique(final_table_up[final_table_up$ID %in% id_up, c("ID", "Description")])

# Save tables of common biological processes
write.table(table_up, file = paste0(output_path_common, '/common_BP_up.txt'), sep = "\t", quote = F, 
            row.names = F, col.names = T) 

####################### UPSET PLOT FOR DOWN-REGULATED GENES ######################
upset_down <- final_table_down %>%
  select(ID, Estres) %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Estres, values_from = Presence, values_fill = list(Presence = 0))

upset_down <- upset_down %>%
  select(ID, sort(colnames(upset_down)[-1]))

stresses = colnames(upset_down)[-1]

upset_down[stresses] = upset_down[stresses] == 1

stresses <- sort(stresses, decreasing = TRUE)

plot_upset_down <- upset(
  upset_down, stresses, name='stresses', width_ratio=0.1, height_ratio=1, sort_sets = FALSE,
  base_annotations=list('Intersection size' = intersection_size(
    text = list(vjust = -0.1, hjust = -0.1, angle = 45, color = "grey"))),
  themes = upset_modify_themes(
    list(
      'main_bar' = theme(text = element_text(size = 20)), 
      'sets' = theme(text = element_text(size = 20)),    
      'intersections_matrix' = theme(text = element_text(size = 15)),
      'intersection_sizes' = theme(text = element_text(size = 20)),  
      'sets_sizes' = theme(text = element_text(size = 20))             
    )
  )
)

# Save upset plots
ggsave(paste0(output_path_common, '/upset_down_plot.png'), plot = plot_upset_down, width = 15, height = 10)

# Identify common down-regulated biological processes in three conditions
id_down <- upset_down %>%
  rowwise() %>%
  filter(sum(c_across(-ID)) == 3) %>%
  pull(ID)

# Extract the unique gene IDs and their descriptions for down-regulated processes
table_down <- unique(final_table_down[final_table_down$ID %in% id_down, c("ID", "Description")])

# Save tables of common biological processes
write.table(table_down, file = paste0(output_path_common, '/common_BP_down.txt'), sep = "\t", quote = F, 
            row.names = F, col.names = T) 

