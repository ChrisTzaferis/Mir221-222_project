# Libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library("ggsci")
library(scales)
library(grid)
library(dplyr)
library(stringr)

kegg_list <- read.delim("KEGG_list.txt")
up_file <- read.delim("UP_TgC_vs_Tg197_KEGG_2019_Mouse_table.txt")
down_file <- read.delim("Down_TgC_vs_Tg197_KEGG_2019_Mouse_table.txt")

kegg_up <- left_join(kegg_list, up_file)
kegg_down <- left_join(kegg_list, down_file)

total_KEGG <- rbind(kegg_up, kegg_down)
total_KEGG$count <- str_count(total_KEGG$Genes, ";") + 1
total_KEGG_filt <- na.omit(total_KEGG)

DF <- total_KEGG_filt
DF$logPval <- as.numeric(-log10(DF$P.value))

ordered_names <- read.delim("KEGG_list.txt") 
level_order <- factor(DF$Term, level = unique(ordered_names$Term))
level_order2 <- factor(DF$Type, level = c("Up", "Down"))

sc <- scale_colour_gradientn(colours = c("deepskyblue", "blue3", "blue4") , na.value = "grey90", oob = squish)

DF2 <- DF

for(i in 1:length(DF2$logPval))
{
  if(DF2$P.value[i] >= 0.05)
  {
    DF2$logPval[i] <- NA
  }
  
  if(DF2$P.value[i] <= 1.577699e-08)
  {
    DF2$logPval[i] <- 4
  }
  
  if(DF2$count[i] < 3)
  {
    DF2$count[i] <- NA
  }
}

png(file = "bubblePlotKegg.png", width = 25, height = 30, units = "cm", res = 300)
p <- ggplot(DF2) + theme_bw() +
  geom_point(aes(x=Type, y=level_order, size=count, color=logPval))+ 
  scale_size(range = c(3, 8.5), breaks = c(4, 6, 8, 10) ) +    
  sc +
  facet_grid(~Type, scales = "free_x") +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
        axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        strip.text.x = element_text(size = 18, color = "Black", face = "bold.italic"),
        strip.background = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid")) +
  labs(x="Type", y="KEGG Pathway", color="-log10(Pval)", size="Number of genes")
p

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("cornflowerblue", "firebrick2")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()
