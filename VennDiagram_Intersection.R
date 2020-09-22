# Venn diagram with VennDiagram package

rm(list=ls()) # Clear previous data

# Load the VennDiagram package
library(VennDiagram)

PathName = setwd(getwd())

# Threshold setting
SCvalueP = 0.8
SCvalueN = -0.8
Pvalue = 0.05

# Color setting
colorsT <- c("#ed652f", "#a332c2", "#9e1b47", "#eb4979", "#cc45ac")

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


####################################################
# Run function: Positive 
SCPos_PV_GeneList_A <- Pvalue_SCPos_Filter(AOri)
SCPos_PV_GeneList_B <- Pvalue_SCPos_Filter(BOri)
SCPos_PV_GeneList_C <- Pvalue_SCPos_Filter(COri)
SCPos_PV_GeneList_D <- Pvalue_SCPos_Filter(DOri)
SCPos_PV_GeneList_E <- Pvalue_SCPos_Filter(EOri)

# Plot the VennDiagram

library("png")

A <- na.omit(SCPos_PV_GeneList_A)
B <- na.omit(SCPos_PV_GeneList_B)
C <- na.omit(SCPos_PV_GeneList_C)
D <- na.omit(SCPos_PV_GeneList_D)
E <- na.omit(SCPos_PV_GeneList_E)


venn.diagram(x = list(A, B, C, D, E) ,
             category.names = c("SARC", "LGG","KIRP","KIRC","LIHC"),
             filename = 'Venn20200922Class1Pos.png',
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
pp <- readPNG("Venn20200922Class1Pos.png")
plot.new() 
rasterImage(pp,0,0,1,1)

# Intersection
Intersect2_AB_Pos <- intersect(SCPos_PV_GeneList_A,SCPos_PV_GeneList_B)
Intersect2_CD_Pos <- intersect(SCPos_PV_GeneList_C,SCPos_PV_GeneList_D)
Intersect4_ABCD_Pos <- intersect(Intersect2_AB_Pos,Intersect2_CD_Pos)
Intersect5_ABCDE_Pos <- intersect(Intersect4_ABCD_Pos,SCPos_PV_GeneList_E)
####################################################

####################################################
# Run function: Negative 
SCNeg_PV_GeneList_A <- Pvalue_SCNeg_Filter(AOri)
SCNeg_PV_GeneList_B <- Pvalue_SCNeg_Filter(BOri)
SCNeg_PV_GeneList_C <- Pvalue_SCNeg_Filter(COri)
SCNeg_PV_GeneList_D <- Pvalue_SCNeg_Filter(DOri)
SCNeg_PV_GeneList_E <- Pvalue_SCNeg_Filter(EOri)

# Plot the VennDiagram

library("png")

A <- na.omit(SCNeg_PV_GeneList_A)
B <- na.omit(SCNeg_PV_GeneList_B)
C <- na.omit(SCNeg_PV_GeneList_C)
D <- na.omit(SCNeg_PV_GeneList_D)
E <- na.omit(SCNeg_PV_GeneList_E)


venn.diagram(x = list(A, B, C, D, E) ,
             category.names = c("SARC", "LGG","KIRP","KIRC","LIHC"),
             filename = 'Venn20200922Class1Neg.png',
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
pp <- readPNG("Venn20200922Class1Neg.png")
plot.new() 
rasterImage(pp,0,0,1,1)

# Intersection
Intersect2_AB_Neg <- intersect(SCNeg_PV_GeneList_A,SCNeg_PV_GeneList_B)
Intersect2_CD_Neg <- intersect(SCNeg_PV_GeneList_C,SCNeg_PV_GeneList_D)
Intersect4_ABCD_Neg <- intersect(Intersect2_AB_Neg,Intersect2_CD_Neg)
Intersect5_ABCDE_Neg <- intersect(Intersect4_ABCD_Neg,SCNeg_PV_GeneList_E)
####################################################

## Export the list
Output_Path = setwd(getwd())
write.table(Intersect5_ABCDE_PosNeg,file=paste0(Output_Path,"/","SCNP",SCvalueP,"_","Intersect5_ABCDE_PosNeg.txt"),row.names = FALSE,col.names = FALSE)
write.table(Intersect5_ABCDE_PosNeg,file=paste0(Output_Path,"/","SCNP",SCvalueP,"_","Intersect5_ABCDE_PosNeg.csv"),row.names = FALSE,col.names = FALSE)

write.table(Intersect5_ABCDE_Pos,file=paste0(Output_Path,"/","SC",SCvalueP,"_","Intersect5_ABCDE_Pos.txt"),row.names = FALSE,col.names = FALSE)
write.table(Intersect5_ABCDE_Pos,file=paste0(Output_Path,"/","SC",SCvalueP,"_","Intersect5_ABCDE_Pos.csv"),row.names = FALSE,col.names = FALSE)

write.table(Intersect5_ABCDE_Neg,file=paste0(Output_Path,"/","SC",SCvalueN,"_","Intersect5_ABCDE_Neg.txt"),row.names = FALSE,col.names = FALSE)
write.table(Intersect5_ABCDE_Neg,file=paste0(Output_Path,"/","SC",SCvalueN,"_","Intersect5_ABCDE_Neg.csv"),row.names = FALSE,col.names = FALSE)