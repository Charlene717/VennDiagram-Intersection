# Venn diagram with VennDiagram package

rm(list=ls()) # Clear previous data

# Load the VennDiagram package
library(VennDiagram)

PathName = setwd(getwd())

# Threshold setting
SCvalueP = 0.8
SCvalueN = -0.8
Pvalue = 0.05
# Files name
AName <- c("/SARC ")
BName <- c("/LGG ")
CName <- c("/KIRP ")
DName <- c("/KIRC ")
EName <- c("/LIHC ")

#################################Function######################################
library("data.table")
# Filter function: Positive and Negative
Pvalue_SCPosNeg_Filter <- function(Cl_Data) {
  Cl_Data_SC <- Cl_Data[Cl_Data$Spearman.s.Correlation >= SCvalueP|Cl_Data$Spearman.s.Correlation <= SCvalueN,]
  Cl_Data_SC_PV <- Cl_Data_SC[Cl_Data_SC$p.Value <= Pvalue,]
  
  Cl_Data_SC_PV_GeneList <- c(as.matrix(Cl_Data_SC_PV[,1]))
  return(Cl_Data_SC_PV_GeneList)
}

## Test
Cl_Data <- SARC
Cl_Data_SCPosNeg_PV_GeneList <- Pvalue_SCPosNeg_Filter(Cl_Data)
## Test

# Filter function: Positive
Pvalue_SCPos_Filter <- function(Cl_Data) {
  Cl_Data_SC <- Cl_Data[Cl_Data$Spearman.s.Correlation >= SCvalueP,]
  Cl_Data_SC_PV <- Cl_Data_SC[Cl_Data_SC$p.Value <= Pvalue,]
  
  Cl_Data_SC_PV_GeneList <- c(as.matrix(Cl_Data_SC_PV[,1]))
  return(Cl_Data_SC_PV_GeneList)
}

# Filter function: Negitive
Pvalue_SCNeg_Filter <- function(Cl_Data) {
  Cl_Data_SC <- Cl_Data[Cl_Data$Spearman.s.Correlation <= SCvalueN,]
  Cl_Data_SC_PV <- Cl_Data_SC[Cl_Data_SC$p.Value <= Pvalue,]
  
  Cl_Data_SC_PV_GeneList <- c(as.matrix(Cl_Data_SC_PV[,1]))
  return(Cl_Data_SC_PV_GeneList)
}

#######################################################################
# Loading the files
APosNeg <- read.csv(paste0(PathName,AName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
BPosNeg <- read.csv(paste0(PathName,BName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
CPosNeg <- read.csv(paste0(PathName,CName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
DPosNeg <- read.csv(paste0(PathName,DName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
EPosNeg <- read.csv(paste0(PathName,EName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)

#######################################################################
# Run function: Positive and Negative
SCPosNeg_PV_GeneList_APosNeg <- Pvalue_SCPosNeg_Filter(APosNeg)
SCPosNeg_PV_GeneList_BPosNeg <- Pvalue_SCPosNeg_Filter(BPosNeg)
SCPosNeg_PV_GeneList_CPosNeg <- Pvalue_SCPosNeg_Filter(CPosNeg)
SCPosNeg_PV_GeneList_DPosNeg <- Pvalue_SCPosNeg_Filter(DPosNeg)
SCPosNeg_PV_GeneList_EPosNeg <- Pvalue_SCPosNeg_Filter(EPosNeg)

# Plot the VennDiagram

library("png")

A <- na.omit(SCPosNeg_PV_GeneList_APosNeg)
B <- na.omit(SCPosNeg_PV_GeneList_BPosNeg)
C <- na.omit(SCPosNeg_PV_GeneList_CPosNeg)
D <- na.omit(SCPosNeg_PV_GeneList_DPosNeg)
E <- na.omit(SCPosNeg_PV_GeneList_EPosNeg)

colorsT <- c("#ed652f", "#c3db0f", "#db750f", "#eb4979", "#cc45ac")
venn.diagram(x = list(A, B, C, D, E) ,
             category.names = c("SARC", "LGG","KIRP","KIRC","LIHC"),
             filename = 'Venn20200922Class1PosNeg.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "Gray",
             fill = colorsT,
             cat.col = colorsT,
             cat.cex = 1,
             margin = 0.15, 
             sub.just =c(1, 1)
)


# Display saved image
options(repr.plot.height=12, repr.plot.width= 12)
library("png")
pp <- readPNG("Venn20200922Class1PosNeg.png")
plot.new() 
rasterImage(pp,0,0,1,1)

# Intersection
Intersect2_AB <- intersect(SCPosNeg_PV_GeneList_APosNeg,SCPosNeg_PV_GeneList_BPosNeg)
Intersect2_CD <- intersect(SCPosNeg_PV_GeneList_CPosNeg,SCPosNeg_PV_GeneList_DPosNeg)
Intersect4_ABCD <- intersect(Intersect2_AB,Intersect2_CD)
Intersect5_ABCDE <- intersect(Intersect4_ABCD,SCPosNeg_PV_GeneList_EPosNeg)
####################################################




