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
AOri <- read.csv(paste0(PathName,AName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
BOri <- read.csv(paste0(PathName,BName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
COri <- read.csv(paste0(PathName,CName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
DOri <- read.csv(paste0(PathName,DName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
EOri <- read.csv(paste0(PathName,EName,"(TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)

#######################################################################
# Run function: Positive and Negative
SCPosNeg_PV_GeneList_A <- Pvalue_SCPosNeg_Filter(AOri)
SCPosNeg_PV_GeneList_B <- Pvalue_SCPosNeg_Filter(BOri)
SCPosNeg_PV_GeneList_C <- Pvalue_SCPosNeg_Filter(COri)
SCPosNeg_PV_GeneList_D <- Pvalue_SCPosNeg_Filter(DOri)
SCPosNeg_PV_GeneList_E <- Pvalue_SCPosNeg_Filter(EOri)

# Plot the VennDiagram

library("png")

A <- na.omit(SCPosNeg_PV_GeneList_A)
B <- na.omit(SCPosNeg_PV_GeneList_B)
C <- na.omit(SCPosNeg_PV_GeneList_C)
D <- na.omit(SCPosNeg_PV_GeneList_D)
E <- na.omit(SCPosNeg_PV_GeneList_E)

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
Intersect2_AB_PosNeg <- intersect(SCPosNeg_PV_GeneList_A,SCPosNeg_PV_GeneList_B)
Intersect2_CD_PosNeg <- intersect(SCPosNeg_PV_GeneList_C,SCPosNeg_PV_GeneList_D)
Intersect4_ABCD_PosNeg <- intersect(Intersect2_AB_PosNeg,Intersect2_CD_PosNeg)
Intersect5_ABCDE_PosNeg <- intersect(Intersect4_ABCD_PosNeg,SCPosNeg_PV_GeneList_E)
####################################################


