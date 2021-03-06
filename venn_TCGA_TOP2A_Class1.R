# Venn diagram with VennDiagram package

rm(list=ls()) # Clear previous data

# Load the VennDiagram package
library(VennDiagram)

FilesName = setwd(getwd())

# Threshold setting
SCvalueP = 0.8
SCvalueN = -0.8
Pvalue = 0.05

#A_SARC
SARC <- read.csv(paste0(FilesName,"/SARC (TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
# Test KIRP_List06 <- subset(KIRP,Cytoband == "1q32.1")
# SARC <- c(as.matrix(SARC))

library("data.table")
SARC_SC08 <- SARC[SARC$Spearman.s.Correlation >= SCvalueP|SARC$Spearman.s.Correlation <= SCvalueN,]
SARC_SC08P005 <- SARC_SC08[SARC_SC08$p.Value <= Pvalue,]

SARC_SC08P005_GeneList <- c(as.matrix(SARC_SC08P005[,1]))

#B_LGG
LGG <- read.csv(paste0(FilesName,"/LGG (TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
LGG_SC08 <- LGG[LGG$Spearman.s.Correlation >= SCvalueP|LGG$Spearman.s.Correlation <= SCvalueN,]
LGG_SC08P005 <- LGG_SC08[LGG_SC08$p.Value <= Pvalue,]

LGG_SC08P005_GeneList <- c(as.matrix(LGG_SC08P005[,1]))

#C_KIRP
KIRP <- read.csv(paste0(FilesName,"/KIRP (TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
KIRP_SC08 <- KIRP[KIRP$Spearman.s.Correlation >= SCvalueP|KIRP$Spearman.s.Correlation <= SCvalueN,]
KIRP_SC08P005 <- KIRP_SC08[KIRP_SC08$p.Value <= Pvalue,]

KIRP_SC08P005_GeneList <- c(as.matrix(KIRP_SC08P005[,1]))

#D_KIRC
KIRC <- read.csv(paste0(FilesName,"/KIRC (TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
KIRC_SC08 <- KIRC[KIRC$Spearman.s.Correlation >= SCvalueP|KIRC$Spearman.s.Correlation <= SCvalueN,]
KIRC_SC08P005 <- KIRC_SC08[KIRC_SC08$p.Value <= Pvalue,]

KIRC_SC08P005_GeneList <- c(as.matrix(KIRC_SC08P005[,1]))

#E_LIHC
LIHC <- read.csv(paste0(FilesName,"/LIHC (TCGA, PanCancer Atlas)_TOP2A.csv"),header = T)
LIHC_SC08 <- LIHC[LIHC$Spearman.s.Correlation >= SCvalueP|LIHC$Spearman.s.Correlation <= SCvalueN,]
LIHC_SC08P005 <- LIHC_SC08[LIHC_SC08$p.Value <= Pvalue,]

LIHC_SC08P005_GeneList <- c(as.matrix(LIHC_SC08P005[,1]))


####################################################

library("png")

A <- na.omit(SARC_SC08P005_GeneList)
B <- na.omit(LGG_SC08P005_GeneList)
C <- na.omit(KIRP_SC08P005_GeneList)
D <- na.omit(KIRC_SC08P005_GeneList)
E <- na.omit(LIHC_SC08P005_GeneList)
#F <- na.omit(LUAD)


## T1
colorsT <- c("#ed652f", "#c3db0f", "#db750f", "#eb4979", "#cc45ac")
venn.diagram(x = list(A, B, C, D, E) ,
             category.names = c("SARC", "LGG","KIRP","KIRC","LIHC"),
             filename = 'Venn20200918Class1SC008V2.png',
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
pp <- readPNG("Venn20200918Class1SC008V2.png")
plot.new() 
rasterImage(pp,0,0,1,1)

####################################################

Intersect2_AB <- intersect(SARC_SC08P005_GeneList,LGG_SC08P005_GeneList)
Intersect2_CD <- intersect(KIRP_SC08P005_GeneList,KIRC_SC08P005_GeneList)
Intersect4_ABCD <- intersect(Intersect2_AB,Intersect2_CD)
Intersect5_ABCDE <- intersect(Intersect4_ABCD,LIHC_SC08P005_GeneList)


