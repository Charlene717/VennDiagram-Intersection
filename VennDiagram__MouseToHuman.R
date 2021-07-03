# Venn diagram with VennDiagram package

rm(list=ls()) # Clear previous data

library(readxl)
library(writexl)

PathName = setwd(getwd())

GSE93601_ALDH2_TtestAll <- read.table(paste0(PathName,"/GSE93601_ALDH2_TtestAll.txt"),  # 資料檔名 
                                      header=T,          # 資料中的第一列，作為欄位名稱
                                      sep="\t")           # 將逗號視為分隔符號來讀取資料

GSE93601_ALDH2_TtestAll_P005 <- GSE93601_ALDH2_TtestAll[GSE93601_ALDH2_TtestAll$Ttest_GeneX_Pvalue < 0.05,]
GSE93601_ALDH2_TtestAll_P005_List <- as.character(GSE93601_ALDH2_TtestAll_P005$ID)

library("tidyverse")
library("readxl")
GO_term1 <- read_excel(paste0(PathName,"/GO_term_summary_20210624_195132.xlsx"))
GO_term2 <- read_excel(paste0(PathName,"/GO_term_summary_20210624_195848.xlsx"))

GO_term_sum <- rbind(GO_term1,GO_term2)
GO_term_sum_list <- GO_term_sum$Symbol
GO_term_sum_list2 <- unique(GO_term_sum_list)
  
# library("dplyr")

# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
genes <- convertMouseGeneList(GO_term_sum_list2)
GO_term_sum_list2_Human <- genes

Intersect2_C12 <- intersect(GSE93601_ALDH2_TtestAll_P005_List,GO_term_sum_list2_Human) 



# Plot the VennDiagram

# Load the VennDiagram package
library(VennDiagram)

colorsT <- c("#e61964", "#f2e22e")
library("png")

A <- na.omit(GSE93601_ALDH2_TtestAll_P005_List)
B <- na.omit(GO_term_sum_list2_Human)


venn.diagram(x = list(A, B) ,
             category.names = c("Acetaldehyde Acc","Immune"),
             filename = 'Venn20201109_V1.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "Gray",
             fill = colorsT,
             #             cat.col = colorsT,
             cat.col = 1,
             cat.cex = 1,
             margin = 0.15, 
             sub.just =c(1, 1)
)


# Display saved image
options(repr.plot.height=12, repr.plot.width= 12)
library("png")
pp <- readPNG("Venn20201109_V1.png")
plot.new() 
rasterImage(pp,0,0,1,1)















##-----------------------------------------------------------------------------
# SumTable <- read.csv(paste0(PathName,"/Gene List of 3 Class.csv"),header = T)
# Class1_List <- SumTable[,2]
# Class2_List <- SumTable[,5]
# Class3_List <- SumTable[,8]
OS  <- read.csv(paste0(PathName,"/PAAD_PValueTableTTT_OSLimit20201118_Xena_15year.csv"))
DFI  <- read.csv(paste0(PathName,"/PAAD_PValueTableTTT_DFILimit20201118_Xena_15year.csv"))
DSS  <- read.csv(paste0(PathName,"/PAAD_PValueTableTTT_DSSLimit20201118_Xena_15year.csv"))
PFI  <- read.csv(paste0(PathName,"/PAAD_PValueTableTTT_PFILimit20201118_Xena_15year.csv"))

OS_List  <- OS[OS$PValue < 0.01,]
OS_List  <- as.character(OS_List[,2])

DFI_List  <- DFI[DFI$PValue < 0.01,]
DFI_List  <- as.character(DFI_List[,2])

DSS_List  <- DSS[DSS$PValue < 0.01,]
DSS_List  <- as.character(DSS_List[,2])

PFI_List  <- PFI[PFI$PValue < 0.01,]
PFI_List  <- as.character(PFI_List[,2])

Intersect2_C12 <- intersect(OS_List,DFI_List)
Intersect2_C34 <- intersect(DSS_List,PFI_List)
Intersect4_C1234 <- intersect(Intersect2_C12,Intersect2_C34)
# Setdiff_C13 <- setdiff(Class1_List,Class3_List)
# Setdiff_C31 <- setdiff(Class3_List,Class1_List)
# 
# Union2_C12 <- union(Class1_List,Class2_List)
# Setdiff_C3UC12 <- setdiff(Class3_List,Union2_C12)
# Setdiff_IC12C3 <- setdiff(Intersect2_C12,Class3_List)
# 
# Union2_C23 <- union(Class2_List,Class3_List)
# Setdiff_C1UC23 <- setdiff(Class1_List,Union2_C23)

# Plot the VennDiagram

# Load the VennDiagram package
library(VennDiagram)

colorsT <- c("#e61964", "#f2e22e", "#29d63a","#2892de")
library("png")

A <- na.omit(OS_List)
B <- na.omit(DFI_List)
C <- na.omit(DSS_List)
D <- na.omit(PFI_List)

venn.diagram(x = list(A, B, C,D) ,
             category.names = c("OS", "DFI","DSS","PFI"),
             filename = 'Venn20201109_V1.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "Gray",
             fill = colorsT,
#             cat.col = colorsT,
             cat.col = 1,
             cat.cex = 1,
             margin = 0.15, 
             sub.just =c(1, 1)
)


# Display saved image
options(repr.plot.height=12, repr.plot.width= 12)
library("png")
pp <- readPNG("Venn20201109_V1.png")
plot.new() 
rasterImage(pp,0,0,1,1)

