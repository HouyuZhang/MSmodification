#=========================================================================================
# This script calculate the intensity of each modification type from Mass-Spectrometry
# result by given a reference table

# Note 1. Please make sure the Mass-Spectrometry returned file was in xlsx format;
# Note 2. Please make sure there aren't irrelevant xlsx files exist in current directory;
# Note 3. You can add desired modification in the "Modification_reference.csv" file.

# Version 1.2 created by Houyu Zhang on 2022/03/18
# Issue report on Hughiez047@gmail.com
# Copyright (c) 2022 __KoziolLab@CIBR__. All rights reserved.
#=========================================================================================
# If a package is installed, it will be loaded. Otherwise the missing package will be installed
# from CRAN and then loaded.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
packages <- c("tidyverse", "readxl", "tools")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressPackageStartupMessages(library(x, character.only = T))
    } else {suppressPackageStartupMessages(library(x, character.only = T))}
  }
)

##' Calculate signal intensity of each modification
##'
##' @param MSFile A excel xlsx records m/z and corresponding signal intensity
##' @param ModificationReferenceFile Specify a output csv file name for standardized-converted lipid signals
##' @examples
##'
##' CalculateModificationIntensity(MSFile = "2022-01-10 Negative-WT1-14 brain regions-03.xlsx",
##'                                ModificationReferenceFile = "Modification_reference.csv")

CalculateModificationIntensity <- function(MSFile = "",
                                           ModificationReferenceFile = ""){

  #Read input files
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

ExcelFiles <- list.files(path = "./", pattern = ".xlsx")
if (length(ExcelFiles) < 1){
  stop("No excel files found in current directory, please check!!!")
  } else {
  cat(length(ExcelFiles),"excel file(s) found in current directory, will calculate one-by-one...\n")
}

#Run each excel file
for (ExcelFile in ExcelFiles){
  CalculateModificationIntensity(MSFile = ExcelFile,
                                 ModificationReferenceFile = "Modification_reference.csv")
}

