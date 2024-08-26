# ==========================================================
#                         01-DMRcaller.R
# ==========================================================
# 
# Description: This script use CX reports of Bismark to extract
# the differentially methylated regions (DMRs) with the bins method.
# Bins method consist in bins where the genome is split and
# all the reads are pooled together
#
# Author: Pascual Villalba y Marta Núñez
# Data: 26/08/24
# Versión: 1.0
#
#
# Notes:
# ----------------------------------------------------------
# 
# Output:
# ----------------------------------------------------------
# 

rm(list = ls())

############################## LIBRARIES #######################################
suppressMessages(library("DMRcaller")) #BiocManager::install("DMRcaller")
suppressMessages(library("argparse"))

################################## FUNCTIONS ###################################

#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = 'DEA.R',
                           description = '
    This program takes the Bismarck CX reports and extract the DMRs',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Project directory path.',
                        required = TRUE)
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'Directory where is the CX report files',
                        required = TRUE)
  required$add_argument('-s', '--summary',
                        type = 'character',
                        help = "Metadata of the study",
                        required = TRUE)
  parser$add_argument('-c', '--cores',
                      default = 1,
                      type = 'double',
                      help = 'cores that the program needs')
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('output', 'input', 'summary','cores')
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
path_out <- args$output
path_CX <- args$input
summary <- args$summary
cores <- args$cores

# Read tables
summary_tab = read.table(summary, sep = "\t", header = F, quote = "\"")
colnames(summary_tab) = c("ID_original", "Time", "Condition", "Replicate")
times = unique(summary_tab$Time)
stresses = unique(summary_tab$Condition)[!grepl("NT", unique(summary_tab$Condition), fixed = T)]

for (t in times) {
  for (s in stresses) {
    cat(paste0("-\tTime: ", t, "; stress: ", s, "\n"))
  
    ## 1-READING CX REPORTS COMING FROM BISMARK
    files = list.files(path_CX, full.names = T)
    files = files[grepl(".CX_report.txt", files, fixed = T)]
    
    Stressed_files = files[grepl(paste0(t, "-", s, "-"), files, fixed = T)]
    NT_files = files[grepl(paste0(t, "-NT-"), files, fixed = T)]
    
    metData_stressed = readBismarkPool(Stressed_files)
    metData_NT = readBismarkPool(NT_files)
    
    metDataList = GRangesList("stressed" = metData_stressed, "NT" = metData_NT)
    
    ## 2-CALLING DMRS
    DMRsBinsCG = computeDMRs(metDataList[["NT"]],
                             metDataList[["stressed"]],
                             regions = NULL,
                             context = "CG",
                             method = "bins",
                             binSize = 50,
                             test = "fisher",
                             pValueThreshold = 0.05,
                             minCytosinesCount = 3,
                             minProportionDifference = 0.1,
                             minGap = 300,
                             minSize = 50,
                             minReadsPerCytosine = 8,
                             cores = cores)
    
    DMRsBinsCHG = computeDMRs(metDataList[["NT"]],
                              metDataList[["stressed"]],
                              regions = NULL,
                              context = "CHG",
                              method = "bins",
                              binSize = 50,
                              test = "fisher",
                              pValueThreshold = 0.05,
                              minCytosinesCount = 3,
                              minProportionDifference = 0.1,
                              minGap = 300,
                              minSize = 50,
                              minReadsPerCytosine = 8,
                              cores = cores)
    
    DMRsBinsCHH = computeDMRs(metDataList[["NT"]], 
                              metDataList[["stressed"]],
                              regions = NULL,
                              context = "CHH",
                              method = "bins",
                              binSize = 50,
                              test = "fisher",
                              pValueThreshold = 0.05,
                              minCytosinesCount = 3,
                              minProportionDifference = 0.1,
                              minGap = 300,
                              minSize = 50,
                              minReadsPerCytosine = 8,
                              cores = cores)
    
    ## 3-MERGING DMRS
    DMRsBinsCGMerged = mergeDMRsIteratively(DMRsBinsCG,
                                            minGap = 200,
                                            respectSigns = TRUE,
                                            metDataList[["NT"]],
                                            metDataList[["stressed"]],
                                            context = "CG",
                                            minProportionDifference = 0.1,
                                            minReadsPerCytosine = 8,
                                            pValueThreshold = 0.05,
                                            test = "fisher",
                                            cores = cores)
    
    DMRsBinsCHGMerged = mergeDMRsIteratively(DMRsBinsCHG,
                                             minGap = 200,
                                             respectSigns = TRUE,
                                             metDataList[["NT"]],
                                             metDataList[["stressed"]],
                                             context = "CHG",
                                             minProportionDifference = 0.1,
                                             minReadsPerCytosine = 8,
                                             pValueThreshold = 0.05,
                                             test = "fisher",
                                             cores = cores)
    
    DMRsBinsCHHMerged = mergeDMRsIteratively(DMRsBinsCHH,
                                             minGap = 200,
                                             respectSigns = TRUE,
                                             metDataList[["NT"]],
                                             metDataList[["stressed"]],
                                             context = "CHH",
                                             minProportionDifference = 0.1,
                                             minReadsPerCytosine = 8,
                                             pValueThreshold = 0.05,
                                             test = "fisher",
                                             cores = cores)
    
    ## 4-WRITING REPORTS
    dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
    
    DMRsBinsCGMerged = as.data.frame(DMRsBinsCGMerged)
    write.table(DMRsBinsCGMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CG.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHGMerged = as.data.frame(DMRsBinsCHGMerged)
    write.table(DMRsBinsCHGMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CHG.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHHMerged = as.data.frame(DMRsBinsCHHMerged)
    write.table(DMRsBinsCHHMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CHH.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    # DMRsBinsCGMerged = DMRsBinsCGMerged[!duplicated(DMRsBinsCGMerged[,1:14]),]
    # write.table(DMRsBinsCGMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CG_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    # DMRsBinsCHGMerged = DMRsBinsCHGMerged[!duplicated(DMRsBinsCHGMerged[,1:14]),]
    # write.table(DMRsBinsCHGMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CHG_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    # DMRsBinsCHHMerged = DMRsBinsCHHMerged[!duplicated(DMRsBinsCHHMerged[,1:14]),]
    # write.table(DMRsBinsCHHMerged, paste0(path_out, t, "-", s, "_DMRs_Bins_CHH_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
  }
}
