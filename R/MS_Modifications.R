#=========================================================================================
# This script contain all available functions for MSI data analysis
# Issue report on zhanghouyu@cibr.ac.cn
# Copyright (c) 2022 __KoziolLab@CIBR__. All rights reserved.
#=========================================================================================

# Loaded installed packages or installed from CRAN and load.
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
packages <- c("tidyverse", "readxl", "tools","MetaboAnalystR",
              "pls","janitor","ggrepel","openxlsx","ggpubr")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = T))
    } else {suppressPackageStartupMessages(library(x, character.only = T))}
  }
)

#=========================================================================================
# Part 1. Functions used for pre-processing MSI rawdata
#=========================================================================================

##' Generate excel tables for store MS-imaging data for each region, each section and each slide
##' @param BrainRegionNames Brain region names selected in the MS-imaging
##' @param SamplesNames Sample names for each section
##' @param SectionNames Number for each section
##' @examples
##'
##'GenrateExcelforMSdata(SamplesNames = c("loxP1","loxP2","loxP3","KO1","KO2","KO3"),
##'                      BrainRegionNames = c("cerebellum","pons&medulla","midbrain","hippocampus","thalamus",
##'                                            "hypothalamus","fornix","caudate putamen","basal forebrain","ventral striatum",
##'                                            "anterior olfactory","olfactory bulbs","cortex","corpus callosum"),
##'                      SectionNames = c("01","02","03"))

GenrateExcelforMSdata <- function(SamplesNames,
                                  BrainRegionNames,
                                  SectionNames){
  for (SamplesName in SamplesNames){
    for (SectionName in SectionNames){
      wb <- createWorkbook(creator = "AUTO")
      for (BrainRegionName in BrainRegionNames){
        addWorksheet(wb = wb, sheetName = paste0(SamplesName,"-",BrainRegionName,"-",SectionName))
      }
      saveWorkbook(wb, file = paste0(Sys.Date()," Negative-",SamplesName,"-14 brain regions-",SectionName,".xlsx"), overwrite = TRUE)
      saveWorkbook(wb, file = paste0(Sys.Date()," Positive-",SamplesName,"-14 brain regions-",SectionName,".xlsx"), overwrite = TRUE)
    }
  }
}

##' Calculate signal intensity of each modification from MS results
##'
##' @param MSFilePath A directory containing excel xlsx files that records m/z and corresponding signal intensity
##' @param ModificationReferenceFile A user defined reference modification list and m/z ranges
##' @examples
##'
##' CalculateModificationIntensity(MSFilePath = "0_Negative rawdata/",
##'                                ModificationReferenceFile = "./Modification_reference_Neg.csv")

CalculateModificationIntensity <- function(MSFilePath,
                                           ModificationReferenceFile){

  MSFiles <- list.files(MSFilePath, pattern = ".xlsx", full.names = T)
  if (length(MSFiles) < 1){
    stop("No excel files found in current directory, please check!!!")
  } else {
    cat(length(MSFiles),"excel file(s) found in current directory, will calculate one-by-one...\n")
  }

  #Run each excel file
  for (MSFile in MSFiles){
    cat("--->>> Processing Excel: [",MSFile,"] ...\n")
    ModificationReference <- read_csv(ModificationReferenceFile, show_col_types = F) %>% as.data.frame()
    MinMR <- min(ModificationReference$mz_lower)
    MaxMR <- max(ModificationReference$mz_upper)
    SheetNames <- excel_sheets(MSFile)

    #Define result matrix format
    ResMat <- matrix(0, nrow = length(unique(ModificationReference$rNs)), ncol = length(SheetNames)) %>% as.data.frame()
    rownames(ResMat) <- unique(ModificationReference$rNs)
    colnames(ResMat) <- SheetNames

    #Process each sheet to get intensity of each modification
    for (SheetName in SheetNames){
      cat("    -> Processing Sheet (",SheetName,") ...\n")
      MSraw <- read_xlsx(MSFile, sheet = SheetName)
      MSraw <- MSraw %>% filter(`m/z` > MinMR & `m/z` < MaxMR)
      MSraw$Base <- ""

      #Process each m/z compound
      if (length(MSraw$`m/z`) > 0){
        for (i in 1:length(MSraw$`m/z`)){
          mz <- MSraw[i,][[1]]
          index <- which(ModificationReference$mz_lower <= mz & ModificationReference$mz_upper >= mz)
          Nucleotide <- ifelse(isTRUE(index > 0), ModificationReference[index,]$rNs, "")
          MSraw[i,]$Base <- Nucleotide
        }
      }
      #Collapse results
      ResTable <- MSraw %>% group_by(Base) %>%
        summarise(Intensity = sum(Intensity)) %>%
        filter(Base != "")

      matchedIndex <- match(ResTable$Base,rownames(ResMat))
      ResMat[matchedIndex,SheetName] <- ResTable$Intensity
    }

    resFileName <- paste0(file_path_sans_ext(MSFile),"_refined.csv")
    cat("--->>> Done Excel: [",MSFile,"] ...\n")
    cat("--->>> Writing refined results to [",resFileName,"] ...\n\n\n")
    write.csv(ResMat, file = resFileName, row.names = T)
  }
}

##' Merged csv files returned by the CalculateModificationIntensity() function
##'
##' @param CSVPath A directory store calculated files generated
##' @param run_RunMetaboAnalystR Whether run MetaboAnalystR analysis for the merged dataset
##' @examples
##'
##' MergeModificationIntensity(CSVPath = "1_Negative refined",
##'                            run_RunMetaboAnalystR = F)

MergeModificationIntensity <- function(CSVPath,
                                       run_RunMetaboAnalystR = F){
  Prefix <- file_path_sans_ext(CSVPath)
  mergedMS <- list.files(path = CSVPath, pattern = ".csv", full.names = T) %>%
    lapply(read_csv) %>%
    bind_cols %>% as.data.frame()

  rownames(mergedMS) <- mergedMS$...1
  mergedMS <- mergedMS %>% select(!starts_with(".."))
  mergedMS <- mergedMS[,order(colnames(mergedMS),decreasing=TRUE)]

  labels <- gsub("-[0-9]+","",colnames(mergedMS))
  mergedMS <- rbind(labels, mergedMS)
  write.csv(file = paste0(Prefix,"-merged.csv"), mergedMS, row.names = T)
  if(run_RunMetaboAnalystR){
    RunMetaboAnalystR(pktablePath = paste0(Prefix,"-merged.csv"), rowNormMet = "NULL")
  }
}

#=========================================================================================
# Part 2. Functions used for normalizing MSI data
#=========================================================================================

##' Normalize signal in each region by dividing total signals in all region in the same section
##'
##' @param MergeModificationIntensityFile The merged modification intensity file generated by MergeModificationIntensity()
##' @param SampleIndex Specify Sample index for group samples
##' @param SectionIndex Specify Section index for group sections
##' @examples
##'
##' NormalizationIntensitySections(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
##'                                SampleIndex = c("loxP1","loxP2","loxP3","KO1","KO2","KO3"),
##'                                SectionIndex = c("-01","-02","-03","-04","-05","-06"))

NormalizationIntensitySections <- function(MergeModificationIntensityFile,
                                           SampleIndex,
                                           SectionIndex){
  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  headers <- read.csv(MergeModificationIntensityFile, header = F, nrows = 1, as.is = T)[-1]
  mergedMS <- read.csv(MergeModificationIntensityFile, header = F, row.names = 1, check.names = F, skip = 2)
  colnames(mergedMS) <- headers

  mergedMSNormaed <- data.frame(matrix(ncol=0, nrow=length(rownames(mergedMS)), dimnames=list(rownames(mergedMS),NULL)))
  #Do normalization for each section of each sample
  for (i in SectionIndex){
    for (j in SampleIndex){
      mergedMS1 <- mergedMS %>% select(contains(i, ignore.case = TRUE)) %>% select(contains(j, ignore.case = TRUE))
      mergedMS1 <- t(apply(mergedMS1,1, function(x) x/sum(x)))
      mergedMSNormaed <- cbind(mergedMSNormaed, mergedMS1)
    }
  }
  labels <- gsub("-[0-9]+","",colnames(mergedMSNormaed))
  mergedMSNormaed <- rbind(labels, mergedMSNormaed)
  write.csv(file = paste0(Prefix,"-NormalizedSections.csv"), mergedMSNormaed, row.names = T)
}

##' Normalize signal in each region by dividing total signals in all region and all sections in the same slide
##'
##' @param MergeModificationIntensityFile The merged modification intensity file generated by MergeModificationIntensity()
##' @param SectionIndex Specify Section index for group sections
##' @examples
##'
##' NormalizationIntensity(MergeModificationIntensityFile = "",
##'                        SectionIndex = c("-01","-02","-03","-04"))

NormalizationIntensitySlides <- function(MergeModificationIntensityFile,
                                         SectionIndex){
  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  headers <- read.csv(MergeModificationIntensityFile, header = F, nrows = 1, as.is = T)[-1]
  mergedMS <- read.csv(MergeModificationIntensityFile, header = F, row.names = 1, check.names = F, skip = 2)
  colnames(mergedMS) <- headers

  mergedMSNormaed <- data.frame(matrix(ncol=0, nrow=length(rownames(mergedMS)), dimnames=list(rownames(mergedMS),NULL)))
  #Do normalization for each section of each sample
  for (i in SectionIndex){
    mergedMS1 <- mergedMS %>% select(contains(i, ignore.case = TRUE))
    mergedMS1 <- t(apply(mergedMS1,1, function(x) x/sum(x)))
    mergedMSNormaed <- cbind(mergedMSNormaed, mergedMS1)
  }
  labels <- gsub("-[0-9]+","",colnames(mergedMSNormaed))
  mergedMSNormaed <- rbind(labels, mergedMSNormaed)
  write.csv(file = paste0(Prefix,"-NormalizedSlides.csv"), mergedMSNormaed, row.names = T)
}

##' Normalize each modification to the total amount of its corresponding nucleotide base
##'
##' @param MergeModificationIntensityFile The merged modification intensity file generated by MergeModificationIntensity()
##' @param ModificationReferenceFile A user defined reference modification list and m/s ranges
##' @examples
##'
##' NormalizationIntensityPercentage(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
##'                                  ModificationReferenceFile = "./Modification_reference_Neg.csv")

NormalizationIntensityPercentage <- function(MergeModificationIntensityFile,
                                             ModificationReferenceFile){

  ModificationReference <- read_csv(ModificationReferenceFile, show_col_types = F) %>% as.data.frame()

  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  headers <- read.csv(MergeModificationIntensityFile, header = F, nrows = 1, as.is = T)[-1]
  mergedMS <- read.csv(MergeModificationIntensityFile, header = F, row.names = 1, check.names = F, skip = 2)
  colnames(mergedMS) <- headers
  mergedMS_Res <- mergedMS
  modificationTypes <- unique(ModificationReference$rNs)
  for (modificationType in modificationTypes){
    BaseGroup <- ModificationReference %>% filter(rNs == modificationType) %>%
      select(Base_type) %>% unique() %>% as.character()
    BaseGroupTypes <- ModificationReference %>% filter(Base_type == BaseGroup) %>% select(rNs) %>% unique()
    mergedMS_Res[modificationType,] <- mergedMS[modificationType,]/colSums(mergedMS[BaseGroupTypes$rNs,])
  }

  labels <- gsub("-[0-9]+","",colnames(mergedMS_Res))
  mergedMS_Res <- rbind(labels, mergedMS_Res)
  write.csv(file = paste0(Prefix,"-NormalizedPercentage.csv"), mergedMS_Res, row.names = T)
}

##' Merge technical replicates based on section suffix
##'
##' @param MergeModificationIntensityFile The merged modification intensity file generated by MergeModificationIntensity()
##' @examples
##'
##' MergeTechnicalReplicates(MergeModificationIntensityFile = "1_Negative refined-merged.csv")

MergeTechnicalReplicates <- function(MergeModificationIntensityFile){

  #Merge technical replicates
  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  mergedMS <- read.csv(MergeModificationIntensityFile, header = T, row.names = 1, check.names = F, skip = 1)

  mergedMS <- t(apply(mergedMS,1, function(x) tapply(x,colnames(mergedMS),mean)))
  labels <- gsub("[1-3]-","-",colnames(mergedMS))
  mergedMS <- rbind(labels, mergedMS)
  write.csv(file = paste0(Prefix,"-mergeTR.csv"), mergedMS, row.names = T)

  #Remove certain sections
  # mergedMSPicked <- mergedMS %>% select(!contains("-01", ignore.case = TRUE))
  # write.csv(file = paste0(Prefix,"-merged_r01.csv"), mergedMSPicked, row.names = T)
}

#=========================================================================================
# Part 3. Functions used for analyzing MSI data
#=========================================================================================

##' Draw modified Volcano Plot using ggplot
##' @param PlotFile The file returned by Volcano.Anal() in MetaboAnalystR
##' @param ThresholdFC Threshold for fold-change
##' @param ThresholdSig Threshold for significance
##' @param TopUpDownShown Select top N up, down significant features for shown
##' @examples
##'
##' volcano_plotting(PlotFile = "1_Negative refined-merged-NormalizedPercentage_volcano.csv",
##'                  ThresholdFC = 1.5,
##'                  ThresholdSig = 0.05,
##'                  TopUpDownShown = c(10,10))

volcano_plotting <- function(PlotFile,
                             ThresholdFC = 1.5,
                             ThresholdSig = 0.05,
                             TopUpDownShown = c(10,10)){
  volcano_pk <- read_csv(PlotFile, show_col_types = FALSE) %>%
    clean_names() %>%
    transform(Group=NA) %>%
    transform(Group=case_when(raw_pval < ThresholdSig & log2_fc > ThresholdFC ~ "Sig.Up",
                              raw_pval < ThresholdSig & log2_fc < -ThresholdFC ~ "Sig.Down",
                              is.na(Group) ~ "Nonsig."))

  tbl <- table(volcano_pk$Group)
  volcano_pk <- volcano_pk %>%
    transform(Group=case_when(Group == "Sig.Up" ~ paste0("Sig.Up [",ifelse(is.na(tbl["Sig.Up"]),0,tbl[["Sig.Up"]]),"]"),
                              Group == "Sig.Down" ~ paste0("Sig.Down [",ifelse(is.na(tbl["Sig.Down"]),0,tbl[["Sig.Down"]]),"]"),
                              Group == "Nonsig." ~ paste0("Nonsig. [",ifelse(is.na(tbl["Nonsig."]),0,tbl[["Nonsig."]]),"]")))

  volcano_pk$Label <- ""

  JudgeNumUp <- ifelse(is.na(tbl["Sig.Up"]),0,tbl[["Sig.Up"]])
  if(JudgeNumUp !=0 ){
    UpShownNum <- ifelse(JudgeNumUp < TopUpDownShown[1], JudgeNumUp, TopUpDownShown[1])
    UpShownItem <- volcano_pk %>% arrange(desc(log10_p)) %>% filter(str_detect(Group, "Sig.Up")) %>% slice(1:UpShownNum)
    volcano_pk$Label[volcano_pk$x1 %in% UpShownItem$x1] <- UpShownItem$x1
  }

  JudgeNumDown <- ifelse(is.na(tbl["Sig.Down"]),0,tbl[["Sig.Down"]])
  if(JudgeNumDown != 0){
    DownShownNum <- ifelse(JudgeNumDown < TopUpDownShown[2], JudgeNumDown, TopUpDownShown[2])
    DownShownItem <- volcano_pk %>% arrange(desc(log10_p)) %>% filter(str_detect(Group, "Sig.Down")) %>% slice(1:DownShownNum)
    volcano_pk$Label[volcano_pk$x1 %in% DownShownItem$x1] <- DownShownItem$x1
  }

  prefix <- file_path_sans_ext(PlotFile)
  pdf(paste0(prefix,"_volcanoPlot_Custome_FC",ThresholdFC,"_Sig",ThresholdSig,".pdf"), height = 8, width = 12)
  p <- ggplot(volcano_pk, aes(x=log2_fc, y=log10_p, color=Group, label=Label)) +
    geom_point(shape=21) +
    theme_bw() + geom_text_repel() +
    scale_color_manual(values=c("grey", "#1F78B4", "#E31A1C")) +
    geom_vline(xintercept=c(-ThresholdFC, ThresholdFC), col="black", linetype="dashed") +
    geom_hline(yintercept=-log10(ThresholdSig), col="black", linetype="dashed") +
    labs(x = "log2(FoldChange)", y = "-log10(P-value)") +
    theme(
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=12, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=12, face="bold"),
      legend.title = element_blank(),
      legend.text = element_text(color="black", size=12, face="bold")
    )
  plot(p)
  dev.off()
}

##' Draw costume scatter plot for pair-wise T.test or multi-group anova result
##' @param PlotFile The file returned by Volcano.Anal()
##' @param ThresholdSig Threshold for significance
##' @param SigIndex select using raw p-value (p_value) or adjusted p-value (fdr)
##' @param TopShown Top significant items for labeling
##' @examples
##'
##' T_Anavo_plotting(PlotFile = "1_Negative refined-merged-NormalizedPercentage_anova_posthoc",
##'                  SigIndex = c("p_value","fdr")[1],
##'                  ThresholdSig = 0.005, TopShown = 10)

T_Anavo_plotting <- function(PlotFile,
                             SigIndex = c("p_value","fdr")[1],
                             ThresholdSig = 0.005,
                             TopShown = 10){
  prefix <- file_path_sans_ext(PlotFile)

  T_Anavo_pk <- read_csv(PlotFile, show_col_types = FALSE) %>% clean_names() %>%
    transform(Group=case_when(eval(parse(text = SigIndex)) < ThresholdSig ~ "Significant",
                              eval(parse(text = SigIndex)) >= ThresholdSig  ~ "Nonsignificant"))

  tbl <- table(T_Anavo_pk$Group)
  T_Anavo_pk <- T_Anavo_pk %>%
    transform(Group=case_when(Group == "Significant" ~ paste0("Significant [",ifelse(is.na(tbl["Significant"]),0,tbl[["Significant"]]),"]"),
                              Group == "Nonsignificant" ~ paste0("Nonsignificant [",ifelse(is.na(tbl["Nonsignificant"]),0,tbl[["Nonsignificant"]]),"]"))) %>%
    arrange(eval(parse(text = SigIndex)))

  T_Anavo_pk$Label <- ""
  JudgeNum <- ifelse(is.na(tbl["Significant"]),0,tbl[["Significant"]])[[1]]

  if(JudgeNum != 0){
    TopShownNum <- ifelse(JudgeNum < TopShown, JudgeNum, TopShown)
    TopShownItem <- T_Anavo_pk %>% filter(str_detect(Group, "Significant")) %>% slice(1:TopShownNum)
    T_Anavo_pk$Label[T_Anavo_pk$x1 %in% TopShownItem$x1] <- TopShownItem$x1
  }

  # T_Anavo_pk$x1 <- factor(T_Anavo_pk$x1, levels = T_Anavo_pk$x1)

  pdf(paste0(prefix,"_TAnavoPlot_Custome_Sig",ThresholdSig,".pdf"), height = 7, width = 10)
  p <- ggplot(T_Anavo_pk, aes(x=x1, y=log10_p, color=Group, label=Label)) +
    geom_point(shape=21) +
    theme_bw() + geom_text_repel() +
    scale_color_manual(values=c("grey", "#1F78B4", "#E31A1C")) +
    geom_hline(yintercept=-log10(ThresholdSig), col="black", linetype="dashed") +
    labs(x = "Features", y = "-log10(P-value)") +
    theme(
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=12, face="bold"),
      legend.title = element_blank(),
      legend.text = element_text(color="black", size=12, face="bold")
    )
  plot(p)
  dev.off()
}


##' Run MetaboAnalystR codes
##' @param pktablePath A standard MetaboAnalystR input file
##' @param rowNormMet method for column normalization, recommend using NULL for data has internal controls
##' @examples
##'
##' RunMetaboAnalystR(pktablePath = "1_Negative refined-merged-NormalizedPercentage.csv",
##'                  rowNormMet = c("SumNorm","NULL")[2])

RunMetaboAnalystR <- function(pktablePath,
                              rowNormMet = c("SumNorm","NULL")[2]){

  prefix <- file_path_sans_ext(pktablePath)
  # Step1. create the mSet Object, specifying that the data to be uploaded
  # is a peak table ("pktable") and that statistical analysis will be performed ("stat").
  mSet <- InitDataObjects(data.type = "pktable", anal.type = "stat", paired = FALSE)
  # read in the filtered peak list
  mSet <- Read.TextData(mSetObj = mSet, filePath = pktablePath, format = "colu", lbl.type = "disc")

  # Step2. perform data processing (filtering/normalization)
  mSet <- SanityCheckData(mSetObj = mSet)
  # Perform data processing - Minimum Value Replacing
  mSet <- ReplaceMin(mSetObj = mSet)
  mSet <- SanityCheckData(mSetObj = mSet)
  mSet <- FilterVariable(mSetObj = mSet, filter = "none", qcFilter = "F", rsd = 20)
  mSet <- PreparePrenormData(mSetObj = mSet)
  mSet <- Normalization(mSetObj = mSet, rowNorm = rowNormMet, transNorm = "NULL", scaleNorm = "NULL", ratio=FALSE, ratioNum=20)
  mSet <- PlotNormSummary(mSetObj = mSet, imgName = paste0(prefix,"_NormFeature_"), format="pdf", dpi = 100, width=NA)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, imgName=paste0(prefix,"_NormSample_"), format="pdf", dpi = 100, width=NA)

  #Step3. Do differential lipid detection
  if(length(unique(mSet$dataSet$prenorm.cls)) == 2){
    cat("Only 2 groups detected, will do t.test...\n")
    mSet <- Ttests.Anal(mSet, nonpar = F, threshp = 0.05, paired = F, equal.var = F, pvalType = "raw", all_results = T)
    # mSet <- PlotTT(mSet, imgName = paste0(prefix,"_Ttest_"), format = "pdf", dpi = 100, width=NA)
    file.rename("fold_change.csv", paste0(prefix,"_fold_change.csv"))
    file.rename("t_test.csv", paste0(prefix,"_t_test.csv"))
    file.rename("t_test_all.csv", paste0(prefix,"_t_test_all.csv"))
    T_Anavo_plotting(PlotFile = paste0(prefix,"_t_test_all.csv"), ThresholdSig = 0.05)

    mSet <- Volcano.Anal(mSet, paired = FALSE, fcthresh = 1, cmpType = 0, nonpar = F, threshp = 1,
                         equal.var = FALSE, pval.type = "raw")
    file.rename("volcano.csv", paste0(prefix,"_volcano.csv"))
    volcano_plotting(PlotFile = paste0(prefix,"_volcano.csv"), ThresholdFC = 1.5, ThresholdSig = 0.05)
    # mSet <- PlotVolcano(mSet, paste0(prefix,"_volcanoPlot_"), plotLbl = 1, format = "pdf", dpi = 100, width=NA)

  } else if (length(unique(mSet$dataSet$prenorm.cls)) > 2){
    cat("More than 2 groups detected, will do Anova test...\n")
    mSet <- ANOVA.Anal(mSetObj = mSet, nonpar = F, thresh = 0.05, post.hoc = "fisher", all_results = FALSE)
    # mSet <- PlotANOVA(mSetObj = mSet, imgName = paste0(prefix,"_anova_"), format = "pdf", dpi = 100, width=NA)
    if(file.exists("anova_posthoc.csv")){
      file.rename("anova_posthoc.csv", paste0(prefix,"_anova_posthoc.csv"))
      T_Anavo_plotting(PlotFile = paste0(prefix,"_anova_posthoc.csv"), ThresholdSig = 0.005)}
  }
  # Step4. Plot overall heatmap view
  mSet <- PlotSubHeatMap(mSet, imgName = paste0(prefix,"_FeatureSignalHeatmap_"), format = "pdf", dpi = 100, width=NA,
                         dataOpt = "norm", scaleOpt = "row", smplDist = "euclidean",clstDist = "ward.D",
                         palette = "bwm", method.nm = "tanova", top.num = 100, viewOpt = "overview",
                         rowV = T, colV = T, border = T, grp.ave = F)
  mSet <- PlotCorrHeatMap(mSet, imgName = paste0(prefix,"_CorrSample_"), format = "pdf", dpi = 100, width=NA,
                          target = "row", cor.method = "pearson", colors = "bwm", viewOpt = "overview", fix.col = T,
                          no.clst = F, corrCutoff = "0")
  file.rename("correlation_table.csv", paste0(prefix,"_CorrSample_table.csv"))
  mSet <- PlotCorrHeatMap(mSet, imgName = paste0(prefix,"_CorrFeature_"), format = "pdf", dpi = 100, width=NA,
                          target = "col", cor.method = "pearson", colors = "bwm", viewOpt = "overview", fix.col = T,
                          no.clst = F, corrCutoff = "0")
  file.rename("correlation_table.csv", paste0(prefix,"_CorrFeature_table.csv"))

  # Step5. perform PCA
  mSet <- PCA.Anal(mSetObj = mSet)
  mSet <- PlotPCAPairSummary(mSetObj = mSet, imgName=paste0(prefix,"_pca_pair_"), format="pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPCAScree(mSetObj = mSet, imgName=paste0(prefix,"_pca_scree_"), format="pdf", dpi = 100, width=NA, scree.num = 5)
  mSet <- PlotPCA2DScore(mSetObj=mSet, imgName=paste0(prefix,"_pca_score2d_"), format="pdf",
                         72, width=NA, pcx=1, pcy=2, reg=0.95, show=1, grey.scale=0)
  mSet <- PlotPCALoading(mSetObj = mSet, imgName=paste0(prefix,"_pca_loading_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCABiplot(mSetObj = mSet, imgName=paste0(prefix,"_pca_biplot_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCA3DScoreImg(mSet, imgName = paste0(prefix,"_pca_score3d_"), "pdf", dpi =100, width=NA, 1,2,3, angl = 40)

  # Step6. perform PLS-DA
  mSet <- PLSR.Anal(mSet, reg=TRUE)
  mSet <- PlotPLSPairSummary(mSet, paste0(prefix,"_pls_pair_"), format = "pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPLS2DScore(mSet, paste0(prefix,"_pls_score2d_"), format = "pdf", dpi = 100, width=NA,
                         inx1 = 1,inx2 = 2,reg = 0.95,show = 1,grey.scale = 0)
  library(pls)
  mSet <- PlotPLS3DScoreImg(mSet, paste0(prefix,"_pls_score3d_"), format = "pdf", dpi = 100, width=NA,1,2,3,40)
  mSet <- PlotPLSLoading(mSet, paste0(prefix,"_pls_loading_"), format = "pdf", dpi = 100, width=NA, 1, 2)
  mSet <- PLSDA.CV(mSet, methodName = "L", compNum = 3, choice = "Q2")
  mSet <- PlotPLS.Classification(mSet, paste0(prefix,"_pls_cv_"), format = "pdf", dpi = 100, width=NA)
  mSet <- PlotPLS.Imp(mSet, paste0(prefix,"_pls_imp_"), format = "pdf", dpi = 100, width=NA,
                      type = "vip", feat.nm = "Comp. 1", feat.num = 20, color.BW = FALSE)

  for (fileName in c("data_orig.qs","preproc.qs","prenorm.qs","row_norm.qs","complete_norm.qs",
                     "pca_loadings.csv","pca_score.csv","plsda_coef.csv","plsda_loadings.csv",
                     "plsda_score.csv","plsda_vip.csv"
  )){
    file.rename(fileName, paste0(prefix,"_",fileName))
  }
}

##' This is the core script for picking certain regions and features for analyses
##'
##' @param MergeModificationIntensityName A directory stores calculated files generated by MergeModificationIntensity()
##' @param PickBrainRegionNames Specify Brain region name(s) for output
##' @param run_RunMetaboAnalystR Whether run MetaboAnalystR analysis in the meanwhile
##' @param FeatureList provide modification types used in the analysis
##' @examples
##'
##' PickBrainRegion(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv",
##'                 PickBrainRegionNames = "cerebellum",
##'                 run_RunMetaboAnalystR = T,
##'                 FeatureList = c("m6Am","ac4C","m22G","hm5CTP","m6dA","m5dC",
##'                                 "ca5dC","m5dCTP","m6dATP","f5dCTP","m5CMP","m6AMP"))

PickBrainRegion <- function(MergeModificationIntensityFile,
                            PickBrainRegionNames,
                            run_RunMetaboAnalystR = T,
                            FeatureList){

  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  mergedMS <- read_csv(MergeModificationIntensityFile) %>% as.data.frame()
  rownames(mergedMS) <- mergedMS$...1
  mergedMS <- mergedMS[,-1]
  # table(gsub("-[0-9]+|[0-9]+-","",colnames(mergedMS))) %>% as.data.frame()

  mergedMSPicked <- mergedMS %>% select(contains(PickBrainRegionNames, ignore.case = TRUE))

  #Defined features for analysis
  if(!is.null(FeatureList)){
    if(!all(FeatureList %in% rownames(mergedMSPicked))){stop("Some features are not in the reference table, please check!!!")}
    FeatureList <- c(rownames(mergedMSPicked)[1],FeatureList)
    mergedMSPicked <- mergedMSPicked[rownames(mergedMSPicked) %in% FeatureList,]
  }
  write.csv(file = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"), mergedMSPicked, row.names = T)

  if(run_RunMetaboAnalystR){
    RunMetaboAnalystR(pktablePath = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"), rowNormMet = "NULL")
  }
}

##' Plot boxplot for each modification within each brain region
##'
##' @param MergeModificationIntensityName A directory stores calculated files generated by MergeModificationIntensity()
##' @param SampleIndex provide the sample prefix to set brain region names
##' @examples
##'
##' BoxModification(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv",
##'                 SampleIndex = c("Neg-mettl3-loxP1-","Neg-mettl3-loxP2-","Neg-mettl3-loxP3-",
##'                                 "Neg-mettl3-KO1-","Neg-mettl3-KO2-","Neg-mettl3-KO3-"))

BoxModification <- function(MergeModificationIntensityFile,
                            SampleIndex){

  Prefix <- file_path_sans_ext(MergeModificationIntensityFile)
  headers <- read.csv(MergeModificationIntensityFile, header = F, nrows = 1, as.is = T)[-1]
  mergedMS <- read.csv(MergeModificationIntensityFile, header = F, row.names = 1, check.names = F, skip = 2)
  colnames(mergedMS) <- headers

  mergedMS <- mergedMS %>% t() %>% as.data.frame()
  mergedMS$SampleName <- rownames(mergedMS)
  mergedMS$BrainRegion <- gsub("-[0-9]+","",mergedMS$SampleName)
  mergedMS$BrainRegion <- gsub(paste0(SampleIndex,collapse = "|"),"",mergedMS$BrainRegion)

  pattern <- gsub("\\[",".*|[",paste0(unique(mergedMS$BrainRegion),collapse = "[0-9]+-"))
  pattern <- paste0("[0-9]+-",pattern,".*")

  mergedMS$Group <- gsub(pattern,"",mergedMS$SampleName)
  mergedMS <- mergedMS %>% clean_names()

  pdf(paste0(Prefix,"_BoxPlot.pdf"), width = 20, height = 4)
  for (i in colnames(mergedMS)[1:c(length(colnames(mergedMS))-3)]){
    Groups <- "group"
    p <- ggplot(mergedMS, aes_string(x = Groups, y = i)) +
      geom_boxplot(aes(fill = group), notch=F, size=0.8) +
      geom_point(size=2, alpha = 0.5) +
      stat_compare_means(method = "t.test", paired = F, comparisons = list(unique(mergedMS$group))) +
      facet_wrap(~brain_region, ncol=14) +
      scale_fill_manual(values = c("#66A61E","#D95F02")) +
      labs(x = "", title = toupper(i)) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=11, face="bold"),
        legend.title = element_text(color="black", size=14, face="bold"),
        legend.text = element_text(color="black", size=12, face="bold"),
        strip.text.x = element_text(size = 10, face="bold")
      )
    plot(p)

  }
  dev.off()
}

