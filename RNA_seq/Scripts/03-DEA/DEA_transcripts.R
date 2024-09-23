################################################################################
##                                                                            
##  DEA.R                                                          
##                                                                            
##  1. Exploratory analysis
##
##  This program takes the tables of absolute counts from a project and
##  performs a Principal Component Analysis (PCA) for each of the time
##  events considered in that project.
##
##  2. Differential expression analysis
##
##  Then, the program performs a differential expression analysis using
##  DESeq2. The absolute counts tables contain a group of control samples and
##  different treatment samples to which they are related. The differential
##  expression analysis is performed considering all possible combinations of
##  control vs treated, so the program returns a result
##  table for each of them. The results table contains all the information
##  provided by the results() function of DESeq2 together with the
##  log2FoldChange and lfcSE from lfcShrink. In addition to the raw data
##  obtained in the analysis, this script also provides tables with those
##  sequences with an adjusted p-value lower than 0.05.
##                                                                            
##                                                                            
##  Author: Marta Núñez Salvador
##  Version: 1.0                                          
##  Date: 18/07/2024                                                          
##                                                                            
###############################################################################

# Load the libraries
suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("argparse"))
suppressMessages(library("tibble"))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("tximport"))

################################## FUNCTIONS ###################################


#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'DEA_transcripts.R',
                           description = '
    This program takes the tables of absolute counts and
    1. Exploratory analysis
    
    This program takes the tables of absolute counts from a project and
    performs a Principal Component Analysis (PCA) for each time of the stress
    events considered in that project.
     
    2. Differential expression analysis
     
    Then, the program performs a differential expression analysis using
    DESeq2. The absolute counts tables contain a group of control samples and
    different treatment samples to which they are related. The differential
    expression analysis is performed considering all possible combinations of
    control vs treated, so the program returns a result table for each of them. 
    The results table contains all the information provided by the results() 
    function of DESeq2 together with the log2FoldChange and lfcSE from lfcShrink. 
    In addition to the raw data obtained in the analysis, this script also provides
    tables with those sequences with an adjusted p-value lower than 0.05.
                           
    3. Volcano plot graphs
                           
    Finally, the program create Volcano Plots for each differential expression
    analysis',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'Salmon result directory',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Differential expression analysis output directory path. If it does not exist, it will be created',
                        required = TRUE)
  required$add_argument('-m', '--metadata',
                        type = 'character',
                        help = 'Additional info directory',
                        required = TRUE)
  parser$add_argument('-a', '--alpha',
                      default = 0.05,
                      type = 'double',
                      help = 'Alpha significance level. Default is 0.05')
  required$add_argument('-s', '--specie',
                        type = 'character',
                        help = "Data's specie name")
  required$add_argument('-p', '--project',
                        type = 'character',
                        help = 'Project name')
  required$add_argument('-g', '--graphs',
                        type = 'character',
                        help = 'Differential expression analysis graphs output directory path. If it does not exist, it will be created')
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('input', 'output', 'metadata','alpha','specie', 'project','graphs')
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
path_in <- args$input
path_metadata <- args$metadata
path_out <- args$output
alpha_value <- args$alpha
specie <- args$specie
project <- args$project
path_graph <- args$graphs

# Create output paths
path_raw_out <- paste(path_out, '01-DEA_raw', specie, project, sep = '/')
path_sig_out <- paste(path_out, '02-DEA_sig', specie, project, sep = '/')
path_out_ea <- paste(path_graph, '01-PCA_graphs', sep = '/')
path_out_vp <- paste(path_graph, '02-Volcano_plots', sep = '/')

# Create directories if they do not exist
dir.create(path_raw_out, recursive = TRUE, showWarnings = FALSE)
dir.create(path_sig_out, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_ea, recursive = TRUE, showWarnings = FALSE)
dir.create(path_out_vp, recursive = TRUE, showWarnings = FALSE)

# Path and tables
metadata <- read.table(paste0(path_metadata,"/metadata.txt"), sep="\t", header=TRUE,row.names = 1, stringsAsFactors = TRUE)
annotation_file <- read.table(paste0(path_metadata,"/DHL92_gene_description_v4.txt"), sep='\t', header= FALSE,quote = "")
colnames(annotation_file) <- c('Gene','Description')

############################# BATCH1 ####################################

metadata_batch1 <- metadata[metadata$Batch == 1,]
files <- file.path(path_in, row.names(metadata_batch1), "quant.sf")
names(files) <- row.names(metadata_batch1)

# This project have three times, each one will be a subproject which will be analised independiently
for (time in unique(metadata_batch1$Time)) {
  
  # Select samples from this subproject
  metadata_subproject <- metadata_batch1[metadata_batch1$Time == time,]
  
  # Obtain samples and files of this group
  files_subproject<- files[rownames(metadata_subproject)]
  
  # The tx2gene file was created using the transcriptome file headers.
  tx2gene <- read.table(paste0(path_metadata,"/tx2gene.txt"), sep = "\t", header = FALSE, stringsAsFactors = TRUE)
  txi <- tximport(files_subproject, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM") 
  
  # Let's construct a DESeqDataSet from the txi `object` and sample information in `sampletable`
  ddsTxi <- DESeqDataSetFromTximport(txi,
                                     colData = metadata_subproject,
                                     design = ~Group)
  # Pre-filtering.
  keep <- rowSums(counts(ddsTxi) > 5) >= 5
  ddsTxi<- ddsTxi[keep,]
  
  # Exploratory analysis and visualization (variance stabilizing transformation)
  vsd_dds <- vst(ddsTxi, blind = FALSE)
  vsd_dds_counts <- assay(vsd_dds)
  
  ## Create and save Principal Component Analysis
  PCA <- plotPCA(vsd_dds, intgroup = c("Group"))
  ggsave(paste0(path_out_ea,"/PCA_time",time,".png"), plot = PCA, width = 8, height = 6, dpi = 300)
  
  # Differential expression analysis
  dds <- DESeq(ddsTxi)
  
  ## Relevel the 'Group' factor to set the specified control group at the given time as the reference level.
  dds$Group <- relevel(dds$Group, ref=paste0("control_",time))
  dds <- DESeq(dds)
  
  ## Obtain results from each contrast
  for(stress in unique(metadata_subproject$Condition)){
    if (stress != "control") {
      ### Extract results for the specified comparison (treatment vs. control at the given time point) 
      res <- results(dds, name=paste0("Group_",stress,"_",time,"_vs_control_",time), alpha = alpha_value)
      
      ### Perform LFC shrinkage to stabilize the estimates, especially for genes with low counts or high variability.
      shrunk <- lfcShrink(dds, coef=paste0("Group_",stress,"_",time,"_vs_control_",time), res=res)
      
      ### Transform Large DESeq Results object to dataframe
      res_tb <- res %>%
        data.frame() %>%
        rownames_to_column(var="seq") %>%
        as_tibble()
      
      ### Add Shrunk values to results table
      res_tb$Shrunkenlog2FoldChange <- shrunk$log2FoldChange
      res_tb$lfcShrunkSE <- shrunk$lfcSE
      
      ### Save table
      write.csv(res_tb,paste0(path_raw_out,"/",stress,"_T",time,"_dea_raw.csv"),row.names = FALSE,quote = FALSE)
      
      #### Extract significant sequences and save it
      sig <- res_tb %>%
        dplyr::filter(padj < alpha_value)
      write.csv(sig,paste0(path_sig_out,"/",stress,"_T",time,"_dea_sig.csv"),row.names = FALSE,quote = FALSE)
      
      ### Create a volcano plot
      res_tb$expression_type <- "No differentially expressed"
      res_tb$expression_type[res_tb$padj < alpha_value & res_tb$Shrunkenlog2FoldChange > 0] <- "UP-regulated"
      res_tb$expression_type[res_tb$padj < alpha_value & res_tb$Shrunkenlog2FoldChange < 0] <- "DOWN-regulated"
      
      #### Calculate the number of sequences of each expression type
      counts <- res_tb %>% 
        group_by(expression_type) %>% 
        summarise(count = n()) %>% 
        ungroup()
      
      #### Create legend labels 
      labels <- counts %>%
        mutate(label = paste0(count, " ", expression_type)) %>%
        pull(label)
      names(labels) <- counts$expression_type
      
      # Calculate limit x
      maxlfc <- max(abs(res_tb$Shrunkenlog2FoldChange))
      
      #### VP
      p <- ggplot(data = res_tb, aes(x = Shrunkenlog2FoldChange, y = -log10(padj), col = expression_type)) + 
        geom_point() + 
        theme_minimal() +
        geom_vline(xintercept = 0, col = "grey") +
        geom_hline(yintercept = -log10(0.05), col = "grey") +
        scale_color_manual(values = c("No differentially expressed" = "snow3", "UP-regulated" = "#FF6F61", "DOWN-regulated" = "#6EC5E9"),
                           labels = labels) +
        labs(title = paste0("Differential Gene Expression in Time ", time, " under ", stress, " stress"), 
             x = "Log2 Fold Change", 
             y = "-Log10 (P-valor ajustado)") +
        scale_x_continuous(limits = c((-maxlfc - 0.5),(maxlfc + 0.5))) +
        theme_bw()
      ggsave(paste0(path_out_vp,"/",stress,"_T",time,".png"), plot = p, width = 8, height = 6, dpi = 300)
    }
  }
}


############################# BATCH2 ####################################

metadata_batch2 <- metadata[metadata$Batch == 2,]
files <- file.path(path_in, row.names(metadata_batch2), "quant.sf")
names(files) <- row.names(metadata_batch2)

# This project have three times, each one will be a subproject which will be analised independiently
for (time in unique(metadata_batch2$Time)) {
  
  # Select samples from this subproject
  metadata_subproject <- metadata_batch2[metadata_batch2$Time == time,]
  
  # Obtain samples and files of this group
  files_subproject<- files[rownames(metadata_subproject)]
  
  # The tx2gene file was created using the transcriptome file headers.
  tx2gene <- read.table(paste0(path_metadata,"/tx2gene.txt"), sep = "\t", header = FALSE, stringsAsFactors = TRUE)
  txi <- tximport(files_subproject, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM") 
  
  # Let's construct a DESeqDataSet from the txi `object` and sample information in `sampletable`
  ddsTxi <- DESeqDataSetFromTximport(txi,
                                     colData = metadata_subproject,
                                     design = ~Group)
  # Pre-filtering.
  keep <- rowSums(counts(ddsTxi) > 5) >= 5
  ddsTxi<- ddsTxi[keep,]
  
  # Exploratory analysis and visualization (variance stabilizing transformation)
  vsd_dds <- vst(ddsTxi, blind = FALSE)
  vsd_dds_counts <- assay(vsd_dds)
  
  ## Create and save Principal Component Analysis
  PCA <- plotPCA(vsd_dds, intgroup = c("Group"))
  ggsave(paste0(path_out_ea,"/PCA_time",time,"batch2.png"), plot = PCA, width = 8, height = 6, dpi = 300)
  
  # Differential expression analysis
  dds <- DESeq(ddsTxi)
  
  ## Relevel the 'Group' factor to set the specified control group at the given time as the reference level.
  dds$Group <- relevel(dds$Group, ref=paste0("control_",time))
  dds <- DESeq(dds)
  
  ## Obtain results from each contrast
  for(stress in unique(metadata_subproject$Condition)){
    if (stress != "control") {
      ### Extract results for the specified comparison (treatment vs. control at the given time point) 
      res <- results(dds, name=paste0("Group_",stress,"_",time,"_vs_control_",time), alpha = alpha_value)
      
      ### Perform LFC shrinkage to stabilize the estimates, especially for genes with low counts or high variability.
      shrunk <- lfcShrink(dds, coef=paste0("Group_",stress,"_",time,"_vs_control_",time), res=res)
      
      ### Transform Large DESeq Results object to dataframe
      res_tb <- res %>%
        data.frame() %>%
        rownames_to_column(var="seq") %>%
        as_tibble()
      
      ### Add Shrunk values to results table
      res_tb$Shrunkenlog2FoldChange <- shrunk$log2FoldChange
      res_tb$lfcShrunkSE <- shrunk$lfcSE
      
      ### Save table
      write.csv(res_tb,paste0(path_raw_out,"/",stress,"_T",time,"_dea_raw.csv"),row.names = FALSE,quote = FALSE)
      
      #### Extract significant sequences and save it
      sig <- res_tb %>%
        dplyr::filter(padj < alpha_value)
      write.csv(sig,paste0(path_sig_out,"/",stress,"_T",time,"_dea_sig.csv"),row.names = FALSE,quote = FALSE)
      
      ### Create a volcano plot
      res_tb$expression_type <- "No differentially expressed"
      res_tb$expression_type[res_tb$padj < alpha_value & res_tb$Shrunkenlog2FoldChange > 0] <- "UP-regulated"
      res_tb$expression_type[res_tb$padj < alpha_value & res_tb$Shrunkenlog2FoldChange < 0] <- "DOWN-regulated"
      
      #### Calculate the number of sequences of each expression type
      counts <- res_tb %>% 
        group_by(expression_type) %>% 
        summarise(count = n()) %>% 
        ungroup()
      
      #### Create legend labels 
      labels <- counts %>%
        mutate(label = paste0(count, " ", expression_type)) %>%
        pull(label)
      names(labels) <- counts$expression_type
      
      # Calculate limit x
      maxlfc <- max(abs(res_tb$Shrunkenlog2FoldChange))
      
      #### VP
      p <- ggplot(data = res_tb, aes(x = Shrunkenlog2FoldChange, y = -log10(padj), col = expression_type)) + 
        geom_point() + 
        theme_minimal() +
        geom_vline(xintercept = 0, col = "grey") +
        geom_hline(yintercept = -log10(0.05), col = "grey") +
        scale_color_manual(values = c("No differentially expressed" = "snow3", "UP-regulated" = "#FF6F61", "DOWN-regulated" = "#6EC5E9"),
                           labels = labels) +
        labs(title = paste0("Differential Gene Expression in Time ", time, " under ", stress, " stress"), 
             x = "Log2 Fold Change", 
             y = "-Log10 (P-valor ajustado)") +
        scale_x_continuous(limits = c((-maxlfc - 0.5),(maxlfc + 0.5))) +
        theme_bw()
      ggsave(paste0(path_out_vp,"/",stress,"_T",time,".png"), plot = p, width = 8, height = 6, dpi = 300)
    }
  }
}

