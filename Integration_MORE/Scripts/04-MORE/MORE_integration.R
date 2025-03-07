################################################################################
##                                                                            ##
## MORE_analysis.R                                                            ##
##                                                                            ##
## 1. Data Processing                                                         ##
##                                                                            ##
## This script processes multiple input files including GO files, metadata,   ##
## methylome associations, and microRNA associations. It integrates these     ##
## datasets to analyze methylome and microRNA interactions in transcriptomic  ##
## data.                                                                      ##
##                                                                            ##
## 2. Analysis Workflow                                                       ##
##                                                                            ##
## The script performs data filtering, normalization, and association         ##
## analyses to identify key interactions between methylation sites,           ##
## microRNAs, and gene expression changes.                                    ##
##                                                                            ##
## 3. Visualization and Export                                                ##
##                                                                            ##
## Results are visualized using various plots and exported as structured      ##
## tables for downstream analysis.                                            ##
##                                                                            ##
## Author: Marta Núñez Salvador                                               ##
## Version: 1.0                                                               ##
## Date: 19/02/2025                                                           ##
##                                                                            ##
################################################################################

suppressMessages(library(argparse))
suppressMessages(library(tibble))
suppressMessages(library(MORE))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.CMelo.eg.db))
################################## FUNCTIONS ###################################

#' Get the command line arguments
#' This function parses the command line arguments entered into the program.
#'
#' @return List with the argument values
get_arguments <- function() {
  parser <- ArgumentParser(prog = 'MORE_analysis.R',
                           description = 'This script analyzes the relationships between methylation, microRNA, and transcriptomic data.')
  
  required <- parser$add_argument_group('required arguments')
  
  required$add_argument('-m', '--metadata', type = 'character',
                        help = 'Path to the metadata file.', required = TRUE)
  required$add_argument('-ma', '--methylome_associations', type = 'character',
                        help = 'Path to the methylome associations file.', required = TRUE)
  required$add_argument('-mi', '--microRNA_associations', type = 'character',
                        help = 'Path to the microRNA associations file.', required = TRUE)
  required$add_argument('-t', '--input_transcripts', type = 'character',
                        help = 'Path to the input transcripts file.', required = TRUE)
  required$add_argument('-im', '--input_microRNA', type = 'character',
                        help = 'Path to the input microRNA file.', required = TRUE)
  required$add_argument('-me', '--input_methylome', type = 'character',
                        help = 'Path to the input methylome file.', required = TRUE)
  required$add_argument('-deg', '--input_DEG', type = 'character',
                        help = 'Path to DEG results.', required = TRUE)
  required$add_argument('-o', '--output', type = 'character',
                        help = 'Path to the output directory.', required = TRUE)
  
  args <- parser$parse_args()
  
  # Check if input files exist
  input_files <- c(args$metadata, args$methylome_associations,
                   args$microRNA_associations, args$input_transcripts,
                   args$input_microRNA, args$input_methylome,args$input_DEG)
  
  for (file in input_files) {
    if (!file.exists(file)) {
      stop(paste('Error: The file', file, 'does not exist.'))
    }
  }
  
  return(args)
}

##################################### MAIN #####################################

# Get program arguments
args <- get_arguments()

# Save the arguments in variables
metadata_path <- args$metadata
methylome_associations_file <- args$methylome_associations
microRNA_associations_file <- args$microRNA_associations
input_transcripts_file <- args$input_transcripts
input_microRNA_file <- args$input_microRNA
input_methylome_file <- args$input_methylome
DEG_path <- args$input_DEG
output_path <- args$output

# Create output directory if it does not exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

############################### MORE ##########################################

# Input tables
transcripts_table <- read.table(input_transcripts_file, sep = "\t", row.names = 1, header = TRUE)
microRNA_table <- read.table(input_microRNA_file, sep = "\t", row.names = 47, header = TRUE)
methylome_upstream_table <- read.table(input_methylome_file,sep = "\t", row.names = 1, header = TRUE)
gene_description <- read.table(paste0(metadata_path,"/DHL92_gene_description_v4.txt"), sep='\t', header= FALSE,quote = "")
colnames(gene_description) <- c('Gene','Description')

# Assosiations table
## microRNA - Genes
association_micro <- read.table(microRNA_associations_file,header = TRUE)
association_micro <- association_micro[,c("gene","microRNA")]

## Methylome - genes
association_methylome_upstream <- read.table(methylome_associations_file,header = FALSE)
colnames(association_methylome_upstream) <- c("Chr", "Start","End","Strand","ID")

## Sorted by chromosome and start position 
association_methylome_upstream <- association_methylome_upstream %>%
  arrange(Chr,Start)


new_rows <- list()
i <- 1  
while (i <= nrow(association_methylome_upstream)) {
  current_row <- association_methylome_upstream[i, ]
  
  # Verify if the next row is < 1000 pb away of the actual
  if (i < nrow(association_methylome_upstream) &&
      association_methylome_upstream$Chr[i] == association_methylome_upstream$Chr[i + 1] &&
      association_methylome_upstream$Start[i + 1] - association_methylome_upstream$Start[i] < 1000) {
    
    # save original row
    new_rows <- append(new_rows, list(current_row))
    
    # Create a copy but with the name of the next gene
    modified_row <- current_row
    modified_row$V5 <- association_methylome_upstream$V5[i + 1]
    
    # Save modify copy
    new_rows <- append(new_rows, list(modified_row))
    
    # skip next row
    i <- i + 1  
  } else {
    # If there is not fusion, save only the actual row
    new_rows <- append(new_rows, list(current_row))
  }
  
  # Come to the next row
  i <- i + 1
}

## Transform list to dataframe
df_association_methylome_upstream<- do.call(rbind, new_rows)

## Create positions names
df_gene_positions <- df_association_methylome_upstream %>%
  mutate(
    Strand = ifelse(Strand == "+", "F", "R"),
    Positions = paste(Chr, Start, End, Strand, sep = "_"),
    gene = ID
  ) %>%
  dplyr::select(gene,Positions)

# Metadata to design MLR
metadata_meth <- read.table(paste0(metadata_path,"/metadata/metadata_methylome.tsv"), sep = "\t",header = TRUE)
metadata_micro <- read.table(paste0(metadata_path,"/metadata/metadata_sRNA.tsv"), sep = "\t",header = TRUE)
metadata_trans <- read.table(paste0(metadata_path,"/metadata/metadata_transcripts.tsv"), sep = "\t",header = TRUE)

# Common samples of the three df
valores_comunes <- Reduce(intersect, list(metadata_micro$MORE, metadata_trans$MORE, metadata_meth$MORE))
metadata_global <- metadata_trans[metadata_trans$MORE %in% valores_comunes,]

# Execute MORE to each stress condition and the respective control
MORE_final <- data.frame()
for (group in unique(metadata_global$Group)){
  if (!grepl("control", group)) {
    print("#########################")
    print(group)
    print("#########################")
    
    # Extract features
    stress_table <-metadata_global[metadata_global$Group == group, c("MORE", "Group")]
    time <- unique(metadata_global[metadata_global$Group == group , "Time"])
    stress <- unique(metadata_global[metadata_global$Group == group , "Condition"])
    control <-metadata_global[metadata_global$Condition == "control" &
                                metadata_global$Time == time, c("MORE", "Group")]
    
    # Create design table
    edesign <- rbind(stress_table, control)
    rownames(edesign) <- NULL
    edesign <- edesign %>%
      column_to_rownames(var = "MORE")
    
    # Filter input tables in function of design table
    transcripts_table_fil <- transcripts_table[, row.names(edesign)]
    microRNA_table_fil <- microRNA_table[, row.names(edesign)]
    methylome_upstream_table_fil <- methylome_upstream_table[, row.names(edesign)]
    methylome_upstream_table_fil <- methylome_upstream_table_fil[row.names(methylome_upstream_table_fil) %in% df_gene_positions$Positions, ]
    
    associations <- list("miRNA-seq" = association_micro,
                         "methylome-upstream" = df_gene_positions)
    
    regulatoryData <- list("miRNA-seq" = microRNA_table_fil,
                           "methylome-upstream" = methylome_upstream_table_fil)
    
    
    file <- paste0(DEG_path,"/",stress,"_T", time,"_dea_sig.csv")
    if (file.exists(file)) {
      
      # Introduce only the genes which are DE in this condition
      DE_genes <- read.table(file,sep = ",", quote = "",header = TRUE)
      genes <- unique(DE_genes$seq)
      transcripts_table_filt_DE <- transcripts_table_fil[rownames(transcripts_table_fil) %in% genes, ]
      if (nrow(transcripts_table_filt_DE) > 1) {
        #####################################
        print(group)
        #####################################
        ## Execute MORE
        SimMLR <- tryCatch({more(
          targetData = transcripts_table_filt_DE,
          associations = associations,
          regulatoryData = regulatoryData,
          condition = edesign,
          scaleType = 'auto',
          varSel = 'EN',
          epsilon = 0.00001,
          alfaEN = NULL,
          interactions = TRUE,
          minVariation = 0,
          correlation = 0.7,
          method  = 'MLR')}, 
                           error = function(e) {
                             message("⚠️ Error en more(): ", e$message)  
                             return(NULL)})
        if (is.null(SimMLR)) {
          print("❌ Error happened...")
          } 
        else {
          print("✅ MORE execute correctly")
          if (length(as.data.frame(SimMLR$GlobalSummary$GoodnessOfFit)) > 1) {
            # Get the results per condition
            MOREregulations <- RegulationPerCondition(SimMLR, filterR2 = 0.9)
          
            MOREregulations_opposite <- MOREregulations[MOREregulations[, 6] < 0 & MOREregulations[, 7] >= 0, ]
            if (nrow(MOREregulations_opposite) != 0) {
              MOREregulations_opposite$Stress <- paste0(stress, "_T", time)
              colnames(MOREregulations_opposite) <- c("Gene", "Regulator", "Omic","Area", "Representative", "Coef.stress", "Coef.control", "Stress")
              
              # Add LFC and padj
              df_DEA_values <- DE_genes[,c("seq","Shrunkenlog2FoldChange","padj")]
              MOREregulations_opposite_des <- merge(merge(MOREregulations_opposite,gene_description, by="Gene"),df_DEA_values, by.x = "Gene",by.y="seq")
            }   
          }
        }
      }
    }
  }
  MORE_final <- rbind(MORE_final, MOREregulations_opposite_des)
} 


MORE_final_sort <- MORE_final %>%
  arrange(Coef.stress)

MORE_final_sort <- unique(MORE_final_sort)

write.table(MORE_final_sort, file = paste0(output_path,"/Regulators_global.tsv"),sep = "\t", col.names = TRUE,row.names = FALSE,quote = FALSE)

# Summarize number of genes with potencial regulators for condition
MORE_final_summary <- MORE_final_sort %>%
  mutate(
    Stress_original = Stress,
    Stress = sub("_.*", "",Stress_original),  
    Time = sub(".*_", "", Stress_original),
    Omic = Omic
  ) %>%
  group_by(Stress, Time, Omic) %>% 
  summarise(
    n_genes = n(),     
    IDs = if_else(Omic == "miRNA-seq", 
                  paste(Gene, Regulator, sep = "-"),
                  Gene) %>% 
      paste(collapse = ", ")
  ) %>%
  ungroup()

write.table(MORE_final_summary, file = paste0(output_path,"/Regulators_summary.tsv"),sep = "\t", col.names = TRUE,row.names = FALSE,quote = FALSE)

