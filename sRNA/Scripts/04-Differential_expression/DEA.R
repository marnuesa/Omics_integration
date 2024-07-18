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
library("DESeq2")
library("ggplot2")
library("dplyr")

# Get programm arguments
args <- get_arguments()

# Save the the arguments in variables
path_table <- args$input
path_metadata <- args$metadata
path_out <- args$output
alpha_value <- args$alpha
specie <- args$specie
project <- args$project
  
# Create output paths
path_raw_out <- paste(path_out, '01-DEA_raw', species, project, sep = '/')
path_sig_out <- paste(path_out, '02-DEA_sig', species, project, sep = '/')
path_out_ea <- paste(path_out, '00-PCA_graphs', sep = '/')
path_out_vp <- paste(path_out, '03-Volcano_plots', sep = '/')

# Create countdata and metadata tables
countdata <- read.csv(path_table, header=TRUE, row.names = "seq",quote = "")
metadata <- read.table(path_metadata, sep=',', header = TRUE, stringsAsFactors = TRUE,row.names = 1)

# This project have three times, each one will be a subproject which will be analised independiently
for (time in unique(metadata$Time)) {
  
  # Select samples from this subproject
  metadata_subproject <- metadata[metadata$Time == time,]
  
  # Filter columns of countdata 
  countdata_subproject <- countdata[,rownames(samplestable_subproject)]
  
  # Create count matrix DESeq input
  dds_matrix <- DESeqDataSetFromMatrix(countData = countdata_subproject,
                                          colData = metadata_subproject,
                                          design = ~ Group)
  
  # Pre-filtering; It filters sequences that are less than 5 times in 5 samples
  keep <- rowSums(countdata_subproject > 5) >= 5
  dds_matrix <- dds_matrix[keep,]
  
  # Exploratory analysis and visualization (variance stabilizing transformation)
  vsd_dds <- vst(dds_matrix, blind = FALSE)
  vsd_dds_counts <- assay(vsd_dds)
  
  ## Create and save Principal Component Analysis
  PCA <- plotPCA(vsd_dds, intgroup = c("Group"))
  ggsave(paste0(path_out_ea,"/PCA_time",time,".png"), plot = PCA, width = 8, height = 6, dpi = 300)
  
  # Differential expression analysis
  dds <- DESeq(dds_matrix)
  
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
      res_tb$expression_type[res_tb$padj < alpha_value] <- "UP-regulated"
      res_tb$expression_type[res_tb$padj < alpha_value] <- "DOWN-regulated"
      
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
      
      #### VP
      p <- ggplot(data = res_tb, aes(x = log2FoldChange, y = -log10(padj), col = expression_type)) + 
        geom_point() + 
        theme_minimal() +
        geom_vline(xintercept = c(-0.6, 0.6), col = "grey") +
        geom_hline(yintercept = -log10(0.05), col = "grey") +
        scale_color_manual(values = c("No differentially expressed" = "snow2", "UP-regulated" = "#FF6F61", "DOWN-regulated" = "#6EC5E9"),
                           labels = labels) +
        labs(title = paste0("Differential Gene Expression in Time ", time, " under ", estres, " stress"), 
             x = "Log2 Fold Change", 
             y = "-Log10 (P-valor ajustado)") +
        scale_x_continuous(limits = c(-10, 10)) +  
        scale_y_continuous(limits = c(0, 10)) +
        theme_bw()
      ggsave(paste0(path_out_vp,"/",estres,"_T",time,".png"), plot = p, width = 8, height = 6, dpi = 300)
    }
  }
}
