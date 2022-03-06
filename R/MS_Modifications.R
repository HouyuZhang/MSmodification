#=========================================================================================
# This script calculate the intensity of each modification type from Mass-Spectrometry
# result by given a reference table

# Note 1. Please make sure the Mass-Spectrometry returned file was in xlsx format;
# Note 2. Please make sure there aren't irrelevant xlsx files exist in current directory;
# Note 3. You can add desired modification in the "Modification_reference.csv" file.

# Version 1.1 created by Houyu Zhang on 2022/03/03
# Issue report on Hughiez047@gmail.com
# Copyright (c) 2022 __KoziolLab@CIBR__. All rights reserved.
#=========================================================================================

# Loaded installed packages or installed from CRAN and load.
packages <- c("tidyverse", "readxl", "tools", "MetaboAnalystR", "pls", "janitor", "ggrepel")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = T))
    } else {suppressPackageStartupMessages(library(x, character.only = T))}
  }
)

BrainRegionNames <- paste0("-",c("cerebellum","pons&medulla","midbrain","hippocampus","thalamus",
                                 "hypothalamus","fornix","caudate putamen","basal forebrain","ventral striatum",
                                 "anterior olfactory","olfactory bulbs","cortex","corpus callosum"))

ModificationsList <- c("methylated A","methylated G","methylated U","methylated C+hm5dC",
                       "I","C","U","T","dI","dC","dA","A+dG",
                       "AMP","CMP","GMP","UMP","TMP","CTP","UTP","ITP","TTP","dCTP","dUTP","dITP","dAMP","dCMP",
                       "m6Am","ac4C","m22G","hm5CTP","m6dA","m5dC","ca5dC","m5dCTP","m6dATP","f5dCTP","m5CMP","m6AMP",
                       "GTP+8-oxo-dGTP","G+8-oxo-dG","m5CTP+hm5dCTP")

##' Draw modified Volcano Plot using ggplot
##' @param PlotFile The file returned by Volcano.Anal() in MetaboAnalystR
##' @param ThresholdFC Threshold for fold-change
##' @param ThresholdSig Threshold for significance
##' @param TopUpDownShown Select top N up, down significant features for shown
##' @examples
##'
##' volcano_plotting(PlotFile = "refinedDF_5depots_abbrev_adi_epi+adi_ibat_volcano.csv",
##'                  ThresholdFC = 1.5,
##'                  ThresholdSig = 0.05,
##'                  TopUpDownShown = c(10,10))

volcano_plotting <- function(PlotFile = "",
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
##' T_Anavo_plotting(PlotFile = "refinedDF_5depots_abbrev_anova_posthoc.csv",
##'                  SigIndex = c("p_value","fdr")[1],
##'                  ThresholdSig = 0.005,TopShown = 10)

T_Anavo_plotting <- function(PlotFile = "",
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
##' RunMetaboAnalystR(pktablePath = "refinedDF_5depots_abbrev_adi_peri+adi_ibat.csv",
##'                 rowNormMet = c("SumNorm","NULL")[1])


RunMetaboAnalystR <- function(pktablePath = "",
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

##' Calculate signal intensity of each modification from MS results
##'
##' @param MSFilePath A directory containing excel xlsx files that records m/z and corresponding signal intensity
##' @param ModificationReferenceFile A user defined reference modification list and m/s ranges
##' @examples
##'
##' CalculateModificationIntensity(MSFile = "2022-01-10 Negative-WT1-14 brain regions-03.xlsx",
##'                                ModificationReferenceFile = "Modification_reference.csv")

CalculateModificationIntensity <- function(MSFilePath = "",
                                           ModificationReferenceFile = ""){

  MSFiles <- list.files(MSFilePath, pattern = ".xlsx")
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
      for (i in 1:length(MSraw$`m/z`)){
        mz <- MSraw[i,][[1]]
        index <- which(ModificationReference$mz_lower <= mz & ModificationReference$mz_upper >= mz)
        Nucleotide <- ifelse(isTRUE(index > 0), ModificationReference[index,]$rNs, "")
        MSraw[i,]$Base <- Nucleotide
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
##' @param PickBrainRegionNames Specify Brain region names for output
##' @param run_RunMetaboAnalystR Whether run MetaboAnalystR analysis in the meanwhile
##' @examples
##'
##' MergeModificationIntensity(CSVPath = "./Negative",
##'                            PickBrainRegionNames = c("hippocampus","cerebellum"))


MergeModificationIntensity <- function(CSVPath = "",
                                       PickBrainRegionNames = c(),
                                       run_RunMetaboAnalystR = T){
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
  # if(run_RunMetaboAnalystR){
  #   RunMetaboAnalystR(pktablePath = paste0(Prefix,"-merged.csv"), rowNormMet = "NULL")
  #   }

  if (!is.null(PickBrainRegionNames)){
    mergedMSPicked <- mergedMS %>% select(contains(PickBrainRegionNames, ignore.case = TRUE))
    write.csv(file = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"), mergedMSPicked, row.names = T)
    if(run_RunMetaboAnalystR){
      RunMetaboAnalystR(pktablePath = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"),
                        rowNormMet = "NULL")
      }
  }
}

##' This is a standalone script for picking certain regions and features
##'
##' @param MergeModificationIntensityName A directory stores calculated files generated by MergeModificationIntensity()
##' @param PickBrainRegionNames Specify Brain region name(s) for output
##' @param run_RunMetaboAnalystR Whether run MetaboAnalystR analysis in the meanwhile
##' @examples
##'
##' PickBrainRegion(MergeModificationIntensityName = "Positive-merged.csv",
##'                 PickBrainRegionNames = "cerebellum",
##'                 run_RunMetaboAnalystR = T,
##'                 FeatureList = c("m6Am","ac4C","m22G","hm5CTP","m6dA","m5dC",
##'                                 "ca5dC","m5dCTP","m6dATP","f5dCTP","m5CMP","m6AMP"))

PickBrainRegion <- function(MergeModificationIntensityName = "",
                            PickBrainRegionNames = c(),
                            run_RunMetaboAnalystR = T,
                            FeatureList = c()){

  Prefix <- file_path_sans_ext(MergeModificationIntensityName)
  mergedMS <- read_csv(MergeModificationIntensityName) %>% as.data.frame()
  rownames(mergedMS) <- mergedMS$...1
  mergedMS <- mergedMS[,-1]

  # mergedMSPicked <- mergedMS %>% select(!contains("-01", ignore.case = TRUE))
  # write.csv(file = paste0(Prefix,"-merged_r01.csv"), mergedMSPicked, row.names = T)

  mergedMSPicked <- mergedMS %>% select(contains(PickBrainRegionNames, ignore.case = TRUE))

  #Defined features for analysis
  if(!is.null(FeatureList)){
    if(!all(FeatureList %in% rownames(mergedMSPicked))){stop("Some features are not in the reference table, please check!!!")}
    FeatureList <- c("1",FeatureList)
    mergedMSPicked <- mergedMSPicked[rownames(mergedMSPicked) %in% FeatureList,]
  }
  write.csv(file = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"), mergedMSPicked, row.names = T)
  if(run_RunMetaboAnalystR){
    RunMetaboAnalystR(pktablePath = paste0(Prefix,paste(PickBrainRegionNames, collapse ="-"),".csv"),
                      rowNormMet = "NULL")
  }
}

