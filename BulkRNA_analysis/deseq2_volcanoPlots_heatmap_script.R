library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)

#--------------------------set directories-------------------------------------- 
setwd(".") 
counts_directory <- "../PROCESSED/" #file containing the counts from FeatureCounts output
sampleFiles <- list.files(path=counts_directory, pattern="*.txt")
sampleTable <- read.delim("Metafile.txt", header=TRUE, sep = "\t")

#---------------------------run deseq2------------------------------------------
ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=counts_directory, design= ~ Condition)

colData(ddsHTseq)$Condition <- factor(colData(ddsHTseq)$Condition, levels=c("WT", "C", "Tg197", "Tg197-C"))
dds <- DESeq(ddsHTseq)

#-------------------------normalization-----------------------------------------
rld <- rlog(dds)
rlogTable <- assay(rld)
write.table(as.data.frame(rlogTable), file="DESeq2_Condition_log-CountTable.txt", row.names = TRUE, col.names = NA ,sep = "\t")

#--------------------------PCA plot--------------------------------------------- 
pcaData <- plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=2.5) +
  theme_bw() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  coord_fixed()
ggsave("condition_PCA.pdf", width = 3.75, height = 2.84, dpi = 600)

#----------------------DESeq2 results-------------------------------------------
control <- c("WT", "WT", "WT")
experimental <- c("C", "Tg197", "Tg197-C")

for(i in 1:length(experimental))
{
  resE_vs_C <- results(dds, contrast=c("Condition", experimental[i], control[i]))
  resE_vs_C <- resE_vs_C[order(resE_vs_C$padj), ]
  file_name <- paste0("DESeq2-Condition_", experimental[i], "_vs_", control[i], ".txt")
  write.table(as.data.frame(resE_vs_C), file=file_name, row.names = TRUE, col.names = NA ,sep = "\t")
}

#-----------------------Volcano plots-------------------------------------------
createVolcano <- function(vecLog2FoldChange, vecPvalue, maxFC, minFC, maxPV, title, fnameOut)
{
  mydata <- data.frame(vecLog2FoldChange, vecPvalue)
  genes <- vector(mode="character", length=length(vecPvalue))
  numUp <- 0
  numDown <- 0
  
  for(i in 1:length(vecPvalue))
  {
    if(vecLog2FoldChange[i] > 1 & vecPvalue[i]<0.05)
    {
      genes[i] <- "Significant up regulated genes"
      numUp <- numUp + 1
    }
    else if(vecLog2FoldChange[i] < -1 & vecPvalue[i]<0.05)
    {
      genes[i] <- "Significant down regulated genes"
      numDown <- numDown + 1
    }
    else
    {
      genes[i] <- "No change"
    }
  }
  
  g <- ggplot(mydata, aes(x=vecLog2FoldChange, y=-log10(vecPvalue))) +
    geom_point(aes(colour = genes), size=0.8)+ #size=5.0
    scale_colour_manual(values = c("Significant up regulated genes"= "red", "Significant down regulated genes"="blue",  "No change"= "black"))+
    xlim(c(minFC, maxFC)) + ylim(c(0, maxPV)) +
    xlab("log2 fold change") + ylab("-log10 Pval") + ggtitle(title) +
    #annotate(geom="text", x = minFC+2, y = maxPV-0.5, label = paste("(", numDown, ")"), color="blue", size=15)+
    #annotate(geom="text", x = maxFC-2, y = maxPV-0.5, label = paste("(", numUp, ")"), color="red", size=15)+
    theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0, vjust = 0),
          axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
          axis.title.y = element_text(face = "bold", color = "black", size = 12),
          axis.title.x = element_text(face = "bold", color = "black", size = 12),
          legend.title = element_text(face = "bold", color = "black", size = 12),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position = "none"
    )
  g
  #g + geom_vline(xintercept = -1, color="red", linetype="dotted", size=2) + geom_vline(xintercept = 1, color="red", linetype="dotted", size=2)
  ggsave(fnameOut, dpi=600, width=2.3, height = 3, units = "in", device = "pdf")
}

wd <- getwd()
DESeq2_files <- list.files(wd)
DESeq2_files_ind <- grep("^DESeq2-Condition", DESeq2_files)
DESeq2_files <- DESeq2_files[DESeq2_files_ind]

for (i in 1:length(DESeq2_files))
{
  file_name <- DESeq2_files[i]
  file_name_interm <- gsub("DESeq2-Condition_", "", file_name)
  file_name_corr <- gsub(".txt", "", file_name_interm)
  file_name_volc <- paste0("volcano_plot_", file_name_corr, ".pdf")
  
  dataF <- read.delim(file_name)
  dataF <- na.omit(dataF)
  
  max_log2_FC <- ceiling((abs(max(dataF$log2FoldChange))))
  min_log2_FC <- ceiling((abs(min(dataF$log2FoldChange))))
  l2FC <- 0
  if(max_log2_FC >= min_log2_FC)
  {
    l2FC <- max_log2_FC
  }
  else
  {
    l2FC <- min_log2_FC
  }
  max_pval <- ceiling((-log10(min(dataF$pvalue))))
  
  createVolcano(dataF$log2FoldChange, dataF$pvalue, l2FC, -l2FC, max_pval,"", file_name_volc)
}


#----------------------filtering part-------------------------------------------
filterDEGs_Up_Down_Pval <- function(file_path, contrast, pval, l2fc_up, l2fc_down)
{
  resE_vs_C <- read.delim(file = file_path, header = T)
  resE_vs_C <- na.omit(resE_vs_C)
  
  deg_resE_vs_C <- resE_vs_C[which(resE_vs_C$pvalue < pval), ]
  upReg_resE_vs_C <- resE_vs_C[which(resE_vs_C$pvalue < pval & resE_vs_C$log2FoldChange > l2fc_up), ]
  downReg_resE_vs_C <- resE_vs_C[which(resE_vs_C$pvalue < pval & resE_vs_C$log2FoldChange < l2fc_down), ]
  
  degs_file_name <- paste0("degs_", contrast, "_pval.txt")
  upReg_file_name <- paste0("up_", contrast, "_pval.txt") 
  downReg_file_name <- paste0("down_", contrast, "_pval.txt")
  
  write.table(deg_resE_vs_C, file=degs_file_name, sep = "\t", col.names = T, row.names = FALSE)
  write.table(upReg_resE_vs_C, file=upReg_file_name, sep = "\t", col.names = T, row.names = FALSE)
  write.table(downReg_resE_vs_C, file=downReg_file_name, sep = "\t", col.names = T, row.names = FALSE)
}

filterDEGs_Up_Down_Pval("DESeq2-Condition_C_vs_WT.txt", "C_vs_WT", 0.05, 1, -1)
filterDEGs_Up_Down_Pval("DESeq2-Condition_Tg197_vs_WT.txt", "Tg197_vs_WT", 0.05, 1, -1)
filterDEGs_Up_Down_Pval("DESeq2-Condition_Tg197-C_vs_WT.txt", "Tg197-C_vs_Wt", 0.05, 1, -1)

#--------------------------Heatmap----------------------------------------------
library(dplyr)
library(ComplexHeatmap)
library(circlize)

degs_1 <- read.delim("degs_C_vs_WT_pval.txt")
degs_2 <- read.delim("degs_Tg197_vs_WT_pval.txt")
degs_3 <- read.delim("degs_Tg197-C_vs_Wt_pval.txt")
all_degs <- rbind(degs_1, degs_2, degs_3)

all_sig <- all_degs[which((all_degs$log2FoldChange > 1 | all_degs$log2FoldChange < -1) & all_degs$pvalue < 0.05), ]

c_vs_wt <- read.delim("DESeq2-Condition_C_vs_WT.txt")
tgc_vs_wt <- read.delim("DESeq2-Condition_Tg197-C_vs_WT.txt")
tg_vs_wt <- read.delim("DESeq2-Condition_Tg197_vs_WT.txt")

all_genes <- as.data.frame(unique(all_sig$X))
colnames(all_genes)[1] <- "X"

final_df_1 <- left_join(all_genes, c_vs_wt)
final_df_1 <- final_df_1[, c(1, 3)]
colnames(final_df_1)[2] <- "C_vs_Wt"

final_df_1 <- left_join(final_df_1, tgc_vs_wt)
final_df_1 <- final_df_1[, c(1, 2, 4)]
colnames(final_df_1)[3] <- "Tg197_C_Vs_Wt"

final_df_1 <- left_join(final_df_1, tg_vs_wt)
final_df_1 <- final_df_1[, c(1, 2, 3, 5)]
colnames(final_df_1)[4] <- "Tg197_Vs_Wt"

final_mat <- as.matrix(final_df_1[, -1])
rownames(final_mat) <- final_df_1$X

clip<-function(x, min=0, max=2) {
  x[x<min]<-min; 
  x[x>max]<-max; 
  x
}

to_plot <- clip(final_mat, -8, 8)
t_to_plot <- t(to_plot)

tiff("heatmap_wt_only.tiff", width = 6.5, height = 3, units = "in", res = 600)
Heatmap(t_to_plot,
        show_row_names = T,
        name = "Log2FC",
        col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
        heatmap_legend_param = list(direction = "vertical"),
        show_column_names = F
)
dev.off()
