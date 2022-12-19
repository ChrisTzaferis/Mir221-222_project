library(Seurat)
library(Scillus)
library(ggplot2)
library(dplyr)

pooled_fibro <- readRDS("Fibro_ATAC.RDS")
pooled_fibro$final_annotation <- factor(pooled_fibro$final_annotation, levels = c("Sublining", "ExpandingIntermediate", "DestructiveLining", "Lining"))
pooled_fibro$seurat_clusters <- pooled_fibro$final_annotation
pooled_fibro$samples2 <- factor(pooled_fibro$samples2, levels = c("wt", "tg"))

#---------------------UMAP visualization----------------------------------------
umap_df <- pooled_fibro@reductions$UMAP@cell.embeddings
umap_df <- as.data.frame(umap_df)
umap_df <- umap_df[, c(1:2)]
umap_df$Cell_id <- rownames(umap_df)
pooled_meta <- pooled_fibro@meta.data
umap_data <- left_join(umap_df, pooled_meta)

mydata=umap_data[, ]
cols = c("aquamarine3", "purple4", "firebrick3", "deepskyblue2")

p <- ggplot(data=mydata ) +
  geom_point(data=mydata[which(mydata$samples2=="tg"), ], aes(x=UMAP_1, y=UMAP_2, fill=final_annotation), size=8, shape=21)+
  scale_fill_manual(values = cols) +
  scale_size()+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        title = element_blank()) +
  labs()
p
ggsave("umap_mouse_tg_ATAC.png", device = "png", width = 11, height = 11, units = "in", dpi = 600)

p <- ggplot(data=mydata ) +
  geom_point(data=mydata[which(mydata$samples2=="wt"), ], aes(x=UMAP_1, y=UMAP_2, fill=final_annotation), size=8, shape=21)+
  scale_fill_manual(values = cols) +
  scale_size()+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        title = element_blank()) +
  labs()
p
ggsave("umap_mouse_wt_ATAC.png", device = "png", width = 11, height = 11, units = "in", dpi = 600)

#-----------------------------Barplot-------------------------------------------
plot_stat(pooled_fibro, plot_type = "prop_fill", group_by = "samples2", pal_setup = cols)+
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 35, angle = 0),
        axis.text.y = element_text(face = "bold", color = "black", size = 35, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 35),
        axis.title.x = element_text(face = "bold", color = "black", size = 35),
        legend.text = element_text(face = "bold", color = "black", size = 30),
        legend.title = element_text(face = "bold", color = "black", size = 30),
        legend.position = "right",
        title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
  labs(x="sample", y="Fraction of cells", fill="Annotation", title = "")
ggsave("barplot_ATAC.png", device = "png", width = 10, height = 15, units = "in", dpi = 600)

#-------------------------Violin plot-------------------------------------------
gene_name <- "Mir221" #Mir222
VlnPlot(pooled_fibro, features = gene_name, split.plot = F, split.by = "samples2", cols = c("grey70", "red"), pt.size = 2, group.by = "final_annotation") +
  ylim(0, 1.5) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 45, angle = 60),
        axis.text.y = element_text(face = "bold", color = "black", size = 45, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 45),
        axis.title.x = element_text(face = "bold", color = "black", size = 45),
        legend.text = element_text(face = "bold", color = "black", size = 40),
        legend.title = element_text(face = "bold", color = "black", size = 40),
        legend.position = "none",
        title = element_text(face = "bold", color = "black", size = 35, angle = 0)) +
  #coord_flip()+
  labs(y="Gene activity score", x="")
ggsave(paste0("violin_", gene_name, "_ATAC.png"), device = "png", width = 12, height = 22, units = "in", dpi = 600)

#-----------------------Feature plot--------------------------------------------
gene_name <- "Mir221"
fp_l <- FeaturePlot(pooled_fibro, features = gene_name, cols = c("grey90", "darkred"), pt.size = 3, order = T, label = F, split.by = "samples2", 
                    label.size = 9, min.cutoff = "q25", max.cutoff = "q95")

ggsave(paste0("FP_", gene_name, "_ATAC.png"), device = "png", width = 21, height = 11, units = "in", dpi = 600)


gene_name <- "Mir222"
fp_l <- FeaturePlot(pooled_fibro, features = gene_name, cols = c("grey90", "navy"), pt.size = 3, order = T, label = F, split.by = "samples2", 
                    label.size = 9, min.cutoff = "q25", max.cutoff = "q95")
ggsave(paste0("FP_", gene_name, "_ATAC.png"), device = "png", width = 21, height = 11, units = "in", dpi = 600)
