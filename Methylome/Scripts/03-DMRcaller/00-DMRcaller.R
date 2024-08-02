################################################################################
#
# Predict DMRs
#
################################################################################

rm(list = ls())

## LIBRARIES
suppressMessages(library(DMRcaller)) #BiocManager::install("DMRcaller")

## VARIABLES
## Create a vector with the arguments.
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("At least 4 arguments must be supplied.", call.=FALSE)
} else {
  WD = args[1]
  CX_path = args[2]
  summary = args[3]
  cores = as.integer(args[4])
}

# WD = "/mnt/doctorado/5-Integracion_omicas_melon/Metiloma/Results/02-DMRcaller/ENDTOEND"
# CX_path = "/mnt/doctorado/5-Integracion_omicas_melon/Metiloma/Results/01-Bismark/ENDTOEND/04-Methylation_extractor"
# summary="/mnt/doctorado/5-Integracion_omicas_melon/Metiloma/Additional_info/Summary_samples/summary_samples.tsv"
# cores = 40

## PIPELINE
summary_tab = read.table(summary, sep = "\t", header = F, quote = "\"")
colnames(summary_tab) = c("ID_original", "Time", "Condition", "Replicate")
times = unique(summary_tab$Time)
stresses = unique(summary_tab$Condition)[!grepl("NT", unique(summary_tab$Condition), fixed = T)]

for (t in times) {
  for (s in stresses) {
    cat(paste0("-\tTime: ", t, "; stress: ", s, "\n"))
  
    ## 1-READING CX REPORTS COMING FROM BISMARK
    files = list.files(CX_path, full.names = T)
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
                             minProportionDifference = 0.15,
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
                              minProportionDifference = 0.15,
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
                              minProportionDifference = 0.15,
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
                                            minProportionDifference = 0.15,
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
                                             minProportionDifference = 0.15,
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
                                             minProportionDifference = 0.15,
                                             minReadsPerCytosine = 8,
                                             pValueThreshold = 0.05,
                                             test = "fisher",
                                             cores = cores)
    
    ## 4-WRITING REPORTS
    if (!dir.exists(paste0(WD, "/DMRs"))){
      dir.create(paste0(WD, "/DMRs"))
    }
    
    DMRsBinsCGMerged = as.data.frame(DMRsBinsCGMerged)
    write.table(DMRsBinsCGMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CG.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHGMerged = as.data.frame(DMRsBinsCHGMerged)
    write.table(DMRsBinsCHGMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CHG.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHHMerged = as.data.frame(DMRsBinsCHHMerged)
    write.table(DMRsBinsCHHMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CHH.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    DMRsBinsCGMerged = DMRsBinsCGMerged[!duplicated(DMRsBinsCGMerged[,1:14]),]
    write.table(DMRsBinsCGMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CG_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHGMerged = DMRsBinsCHGMerged[!duplicated(DMRsBinsCHGMerged[,1:14]),]
    write.table(DMRsBinsCHGMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CHG_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    DMRsBinsCHHMerged = DMRsBinsCHHMerged[!duplicated(DMRsBinsCHHMerged[,1:14]),]
    write.table(DMRsBinsCHHMerged, paste0(WD, "/DMRs/", t, "-", s, "_DMRs_Bins_CHH_duplicates_removed.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
}

