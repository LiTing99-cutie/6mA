###### 1. generate rawcounts and fpkm matrix ######
if(F){

  setwd("/home/user/data2/lit/project/6mA/RNA-seq/")
  
  rm(list=ls())
  
  library(tibble) # unite
  library(stringr) # str_split_fixed
  library(tidyverse)
  source("/home/user/data2/lit/project/1KMG/RNA_seq_cross_4_species/myUtils.R")
  
  counts <- read.table("./output/featureCounts/featureCounts.txt",header=T,skip=1,sep="\t")
  counts <- counts[,c(1,7:ncol(counts))]
  colnames(counts) <- c("WB_gene_id",
                        colnames(counts)[2:ncol(counts)] %>% str_split_fixed(.,"\\.",n=Inf) %>% subset(.,select=c(12)))
  counts %<>% column_to_rownames("WB_gene_id")
  counts.1 <- counts
  
  counts <- read.table("./output/featureCounts/add/featureCounts.txt",header=T,skip=1,sep="\t")
  counts <- counts[,c(1,7:ncol(counts))]
  colnames(counts) <- c("WB_gene_id",
                        colnames(counts)[2:ncol(counts)] %>% str_split_fixed(.,"\\.",n=Inf) %>% subset(.,select=c(12)))
  counts %<>% column_to_rownames("WB_gene_id")
  counts.2 <- counts
  
  cbind(counts.1,counts.2) -> counts
  
  write.table(counts,file = "./output/featureCounts/rawcounts.csv",sep=",")
  
  
  libSize <- read.table("output/rpkm/mappedReadCounts.txt",header = F,sep = "\t",stringsAsFactors = F)
  colnames(libSize) <- c("sample","libSize")
  all(libSize$sample == colnames(counts))
  
  
  
  
  genLen <- read.table("./output/rpkm/gene_length/ce11.ensembl.105.geneLen.txt",header=T,row.names =1) 
  
  t <- match(rownames(genLen),rownames(counts))
  counts.o <- counts[t,]
  
  
  fpkm_1=edgeR_rpkm(counts.o,libSize,genLen,geneLenCol = "merged",
                    group=c(rep("N2_OP50",2),rep("N2_PA14",2),
                            rep("KO_OP50",2),rep("KO_PA14",2),
                            rep("M20_OP50",2),rep("M20_PA14",2)),normLibSize = T,normLog2 = F)
  
  write.table(fpkm_1,file = "./output/featureCounts/fpkm_1.csv",sep=",")
  
  
  genLen <- read.table("./output/featureCounts/featureCounts.txt",header=T,skip=1,sep="\t") %>% .$Length %>% data.frame() 
  colnames(genLen) <- "featureCounts"
  rownames(genLen) <-  rownames(counts)
  
  fpkm=edgeR_rpkm(counts,libSize,genLen,geneLenCol = "featureCounts",
                  group=c(rep("N2_OP50",2),rep("N2_PA14",2),
                          rep("KO_OP50",2),rep("KO_PA14",2),
                          rep("M20_OP50",2),rep("M20_PA14",2)),normLibSize = T,normLog2 = F)
  
  write.table(fpkm,file = "./output/featureCounts/fpkm.csv",sep=",")
  

}


###### 2. DEG ######


library(DESeq2)

mycounts <- read.table("./output/featureCounts/rawcounts.csv", header = T, row.names = 1,sep = ',')

head(mycounts)


idTrans <- read.table("/home/user/data/lit/database/public/annotation/gene_and_gene_predictions/ce11.wormbaseID2geneName.txt", 
                      header = T, sep = ',')

t <- match(rownames(mycounts),idTrans$WormBase.Gene.ID)

rownames(mycounts) <- idTrans[t,]$Gene.name 

condition <- factor(c(rep("N2_OP50",2),rep("N2_PA14",2),
                      rep("KO_OP50",2),rep("KO_PA14",2),
                      rep("M20_OP50",2),rep("M20_PA14",2)), levels = c("N2_OP50","N2_PA14","KO_OP50",
                                                                       "KO_PA14","M20_OP50","M20_PA14"))

condition

colData <- data.frame(row.names = colnames(mycounts), condition)

colData

dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)

dds <- DESeq(dds,parallel = TRUE)

dds

res = results(dds, contrast=c("condition", "N2_PA14", "N2_OP50"))
head(res)
diff_gene_PA_N2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_gene_PA_N2_unfil <- res
dim(diff_gene_PA_N2)
head(diff_gene_PA_N2)

res = results(dds, contrast=c("condition", "KO_PA14", "KO_OP50"))
head(res)
diff_gene_PA_KO <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_gene_PA_KO_unfil <- res
dim(diff_gene_PA_KO)
head(diff_gene_PA_KO)

res = results(dds, contrast=c("condition", "M20_PA14", "M20_OP50"))
head(res)
diff_gene_PA_M20 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_gene_PA_M20_unfil <- res
dim(diff_gene_PA_M20)
head(diff_gene_PA_M20)

write.table(data.frame(diff_gene_PA_N2_unfil),file = "./output/DEG/in_house.DEG.PA.in.WT.csv",sep = ',',row.names = T,col.names = T)
write.table(data.frame(diff_gene_PA_KO_unfil),file = "./output/DEG/in_house.DEG.PA.in.KO.csv",sep = ',',row.names = T,col.names = T)
write.table(data.frame(diff_gene_PA_M20_unfil),file = "./output/DEG/in_house.DEG.PA.in.M20.csv",sep = ',',row.names = T,col.names = T)


###### 3. intersect ######

WT <- read.table("./output/DEG/ref_LY/DEG.PA.in.WT.csv",sep = ',',header = T,row.names = 1)
KO <- read.table("./output/DEG/ref_LY/DEG.PA.in.KO.csv",sep = ',',header = T,row.names = 1)
M20 <- read.table("./output/DEG/ref_LY/DEG.PA.in.M20.csv",sep = ',',header = T,row.names = 1)

xx <- function(d1,d2){
  dd1 <- length(rownames(d1))
  dd2 <- length(rownames(d2))
  xxx <- length(intersect(rownames(d1),rownames(d2)))
  per <- length(intersect(rownames(d1),rownames(d2)))/min(dd1,dd2)
  return(c(dd1,dd2,xxx,per))
}

xx(WT,diff_gene_PA_N2)
xx(KO,diff_gene_PA_KO)
xx(M20,diff_gene_PA_M20)

install.packages("VennDiagram")
library(VennDiagram)

xx <- function(d1,d2,file){
  dd1 <- length(rownames(d1))
  dd2 <- length(rownames(d2))
  xxx <- length(intersect(rownames(d1),rownames(d2)))
  venn.plot <- draw.pairwise.venn(
    area1 = dd1,  #区域1的数
    area2 = dd2,   #区域2的数
    cross.area = xxx,  #交叉数
    category = c("LY","in_house"),#分类名称
    fill = c("blue", "pink"), #区域填充颜色
    cat.dist =0,
    rotation.degree = 0
  )
  grid.newpage()
  png(file)
  grid.draw(venn.plot) #显示图形
  dev.off()
}

xx(WT,diff_gene_PA_N2,"WT.png")
xx(KO,diff_gene_PA_KO,"KO.png")
xx(M20,diff_gene_PA_M20,"M20.png")

       