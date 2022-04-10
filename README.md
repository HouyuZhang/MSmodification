# :hammer_and_wrench:MSmodification: Mass spectrometry-based measurement of nucleotide modifications in brain system

**:writing_hand:Author**: Houyu Zhang and Yan Peng   **:email:Email**: zhanghouyu@cibr.ac.cn

Copyright (c) 2022 KoziolLab@CIBR. All rights reserved.

### Introduction

[Magdalena J Koziol Lab](https://dnalaboratory.org/) is studying and aiming to discover further novel DNA and RNA modifications and investigate their role in the brain using mass spectrometry, MS-imaging (MSI) and nanopore sequencing, etc. 

The MSmodification package is developed to analyze the MSI data.

### Quick start

1. **Installation and load the package**

```R
#MSmodification package can be installed from github using:
devtools::install_github("HouyuZhang/MSmodification")
#Load the installed package
library(MSmodification)
```

2. **MSI raw data preprocessing**

- Calculate the modification intensity based on the reference table, you can find reference tables [here](https://github.com/HouyuZhang/MSmodification/blob/master/CustomeScripts/MSmodification_referenceList.xlsx).

```R
CalculateModificationIntensity(MSFilePath = "0_Negative rawdata/",
                               ModificationReferenceFile = "./Modification_reference_Neg.csv")
```

- Merge each samples into a consensus table:

```R
MergeModificationIntensity(CSVPath = "1_Negative refined", 
                           run_RunMetaboAnalystR = F)
```

- Normalize the intensity using three methods. 

```R
#Normalization across each slide
NormalizationIntensitySlides(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                             SectionIndex = c("-01","-02","-03","-04","-05","-06"))
#Normalization using the percentage
NormalizationIntensityPercentage(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                                 ModificationReferenceFile = "./Modification_reference_Neg.csv")
#Normalization across each section
NormalizationIntensitySections(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                               SampleIndex = c("loxP1","loxP2","loxP3","KO1","KO2","KO3"),
                               SectionIndex = c("-01","-02","-03","-04","-05","-06"))
```

- ​	This package provide four types of data, there are some suggestions on selecting data types for downstream applications: 

  - `Rawdata`: The rawdata is highly meaningful and straightforward because it is calculated by referring to standards, however you need to carefully compare the values across different slides.

  - `Sections-based`: Not recommended, unless you only compare the results within the same section.
  - `Slide-based`: Recommended when comparing data across slides and mapping the intensity to MS images. The relative values in the MS images (blue and red color in MassImager) is equivalent to the slide-based normalized values.
  - `Percentage-based`: Recommend when you want to check whether the modification ratio is different, because higher total nucleoside concentration can results higher modified nucleoside without mechanical changes.

- Analyze modification dynamics based on selected brain regions and modification types

```R
PickBrainRegion(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv",
                PickBrainRegionNames = "cerebellum",
                run_RunMetaboAnalystR = T,
                FeatureList = c("m6Am","ac4C","m22G","hm5CTP","m6dA","m5dC",
                                "ca5dC","m5dCTP","m6dATP","f5dCTP","m5CMP","m6AMP"))
```

​	4. Analyze each modification in each brain region region

```R
BoxModification(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv",
                SampleIndex = c("mettl3-loxP1-","mettl3-loxP2-","mettl3-loxP3-",
                                "mettl3-KO1-","mettl3-KO2-","mettl3-KO3-"))
```

### How to generate the MSI rawdata?

- The whole mouse brain is dissected and frozen at -80 degree, then brain was sliced into 20μm `sections` covering the `sides`

- Slides was subjected to the MSI machine for scanning with a full channel

- MSI results can be visualized by MassImager, brain regions of interest can be segmented to extract the m/z information

- Using the built-in function to generate excel with sheet to store MSI segmentation results. I will always suggest using this function since misspells affect downstream analyses

  ```R
  GenrateExcelforMSdata(SamplesNames = c("loxP1","loxP2","loxP3","KO1","KO2","KO3"),
                        BrainRegionNames = c("cerebellum","pons&medulla","midbrain",
                                             "hippocampus","thalamus","hypothalamus","fornix",
                                             "caudate putamen","basal forebrain",
                                             "ventral striatum","anterior olfactory",
                                             "olfactory bulbs","cortex","corpus callosum"),
                        SectionNames = c("01","02","03"))
  ```

