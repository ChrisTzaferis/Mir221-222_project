######################################
############### Tracks ###############
######################################

## Load object ##
fibro6a_Co_NEW <- loadArchRProject("/media/raid/users/dimitris/data/scATAC-seq/Archr_final/fibro6a_Co_NEW_new/")

## Create new object ##
fibro6a_Co_NEW_mir_tg<-subsetArchRProject(ArchRProj = fibro6a_Co_NEW,cells = fibro6a_Co_NEW[fibro6a_Co_NEW$samples2=="tg",]$cellNames, outputDirectory = "fibro6a_Co_NEW_mir_tg_new",dropCells=FALSE, force=T)
fibro6a_Co_NEW_mir_wt<-subsetArchRProject(ArchRProj = fibro6a_Co_NEW,cells = fibro6a_Co_NEW[fibro6a_Co_NEW$samples2=="wt",]$cellNames, outputDirectory = "fibro6a_Co_NEW_mir_wt_new",dropCells=FALSE, force=T)

## Create new peakSet - add p2g links ##
fibro6a_Co_NEW_mir_tg <- addGroupCoverages(ArchRProj = fibro6a_Co_NEW_mir_tg, groupBy = "Clusters4", force=T)
fibro6a_Co_NEW_mir_tg <- addReproduciblePeakSet(ArchRProj = fibro6a_Co_NEW_mir_tg, groupBy = "Clusters4", pathToMacs2 = findMacs2(), force=T, reproducibility = "1")
fibro6a_Co_NEW_mir_tg <- addPeakMatrix(fibro6a_Co_NEW_mir_tg,force=T)

fibro6a_Co_NEW_mir_tg <- addPeak2GeneLinks(
  ArchRProj = fibro6a_Co_NEW_mir_tg,
  reducedDims = "IterativeLSI"
  , useMatrix = "GeneScoreMatrix")

fibro6a_Co_NEW_mir_wt <- addGroupCoverages(ArchRProj = fibro6a_Co_NEW_mir_wt, groupBy = "Clusters4", force=T)
fibro6a_Co_NEW_mir_wt <- addReproduciblePeakSet(ArchRProj = fibro6a_Co_NEW_mir_wt, groupBy = "Clusters4", pathToMacs2 = findMacs2(), force=T, reproducibility = "1")
fibro6a_Co_NEW_mir_wt <- addPeakMatrix(fibro6a_Co_NEW_mir_wt,force=T)

fibro6a_Co_NEW_mir_wt <- addPeak2GeneLinks(
  ArchRProj = fibro6a_Co_NEW_mir_wt,
  reducedDims = "IterativeLSI"
  , useMatrix = "GeneScoreMatrix")

## Peaks with motifs ##

fibro6a_Co_NEW_mir_tg <- addMotifAnnotations(ArchRProj = fibro6a_Co_NEW_mir_tg, motifSet = "cisbp", name = "Motif",force = TRUE)
fibro6a_Co_NEW_mir_tg <- addBgdPeaks(fibro6a_Co_NEW_mir_tg)
fibro6a_Co_NEW_mir_tg <- addDeviationsMatrix(
  ArchRProj = fibro6a_Co_NEW_mir_tg, 
  peakAnnotation = "Motif",
  force = TRUE
)

fibro6a_Co_NEW_mir_wt <- addMotifAnnotations(ArchRProj = fibro6a_Co_NEW_mir_wt, motifSet = "cisbp", name = "Motif",force = TRUE)
fibro6a_Co_NEW_mir_wt <- addBgdPeaks(fibro6a_Co_NEW_mir_wt)
fibro6a_Co_NEW_mir_wt <- addDeviationsMatrix(
  ArchRProj = fibro6a_Co_NEW_mir_wt, 
  peakAnnotation = "Motif",
  force = TRUE
)

lining_annotation1<-read.delim("/media/raid/users/dimitris/data/scATAC-seq/Archr_final/S4a_subclustering/lining_annotation1.txt",header=F)
lining_annotation1_S4a0_cells<-lining_annotation1[which(lining_annotation1$V2 %in% "S4a_0"),]$V1
lining_annotation1_S4a1_cells<-lining_annotation1[which(lining_annotation1$V2 %in% "S4a_1"),]$V1

super_clusters <- fibro6a_Co_NEW_mir_tg$predictedGroup
super_clusters[which(super_clusters %in% c("S1","S2.a","S2.b","S2.c","S3","S5"))]<-"homeostatic"
super_clusters[which(super_clusters %in% c("S2.d","S4.b"))]<-"intermediate"
super_clusters[which(super_clusters %in% c("S4.a"))]<-"lining"
names(super_clusters)<-rownames(fibro6a_Co_NEW_mir_tg)
super_clusters[which(names(super_clusters) %in% lining_annotation1_S4a0_cells)]<-"S4.a_0"
super_clusters[which(names(super_clusters) %in% lining_annotation1_S4a1_cells)]<-"S4.a_1"
super_clusters2<-paste0(super_clusters,"_",fibro6a_Co_NEW_mir_tg$samples2)
fibro6a_Co_NEW_mir_tg <- addCellColData(ArchRProj = fibro6a_Co_NEW_mir_tg, data = super_clusters,cells = fibro6a_Co_NEW_mir_tg$cellNames, name = "super_clusters",force=T)
fibro6a_Co_NEW_mir_tg <- addCellColData(ArchRProj = fibro6a_Co_NEW_mir_tg, data = super_clusters2,cells = fibro6a_Co_NEW_mir_tg$cellNames, name = "super_clusters2",force=T)
saveArchRProject(fibro6a_Co_NEW_mir_tg, outputDirectory = "fibro6a_Co_NEW_mir_tg_new", load = FALSE)

super_clusters <- fibro6a_Co_NEW_mir_wt$predictedGroup
super_clusters[which(super_clusters %in% c("S1","S2.a","S2.b","S2.c","S3","S5"))]<-"homeostatic"
super_clusters[which(super_clusters %in% c("S2.d","S4.b"))]<-"intermediate"
super_clusters[which(super_clusters %in% c("S4.a"))]<-"lining"
names(super_clusters)<-rownames(fibro6a_Co_NEW_mir_wt)
super_clusters[which(names(super_clusters) %in% lining_annotation1_S4a0_cells)]<-"S4.a_0"
super_clusters[which(names(super_clusters) %in% lining_annotation1_S4a1_cells)]<-"S4.a_1"
super_clusters2<-paste0(super_clusters,"_",fibro6a_Co_NEW_mir_wt$samples2)
fibro6a_Co_NEW_mir_wt <- addCellColData(ArchRProj = fibro6a_Co_NEW_mir_wt, data = super_clusters,cells = fibro6a_Co_NEW_mir_wt$cellNames, name = "super_clusters",force=T)
fibro6a_Co_NEW_mir_wt <- addCellColData(ArchRProj = fibro6a_Co_NEW_mir_wt, data = super_clusters2,cells = fibro6a_Co_NEW_mir_wt$cellNames, name = "super_clusters2",force=T)
saveArchRProject(fibro6a_Co_NEW_mir_wt, outputDirectory = "fibro6a_Co_NEW_mir_tg_new", load = FALSE)

motifMatrix_fibro6a_Co_NEW_mir_tg<-SummarizedExperiment::assay(getMatches(fibro6a_Co_NEW_mir_tg))
rownames(motifMatrix_fibro6a_Co_NEW_mir_tg) <- paste0(seqnames(getMatches(fibro6a_Co_NEW_mir_tg)), ":", start(getMatches(fibro6a_Co_NEW_mir_tg)), "-", end(getMatches(fibro6a_Co_NEW_mir_tg)))
motifMatrix_fibro6a_Co_NEW_mir_tg <- as(motifMatrix_fibro6a_Co_NEW_mir_tg, "dgCMatrix")

motifMatrix_fibro6a_Co_NEW_mir_tg<-SummarizedExperiment::assay(getMatches(fibro6a_Co_NEW_mir_tg))
rownames(motifMatrix_fibro6a_Co_NEW_mir_tg) <- paste0(seqnames(getMatches(fibro6a_Co_NEW_mir_tg)), ":", start(getMatches(fibro6a_Co_NEW_mir_tg)), "-", end(getMatches(fibro6a_Co_NEW_mir_tg)))
motifMatrix_fibro6a_Co_NEW_mir_tg <- as(motifMatrix_fibro6a_Co_NEW_mir_tg, "dgCMatrix")

Relb_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(856)]
Relb_peaks_tg<-names(Relb_peaks_tg[Relb_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Relb_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Relb_peaks_tg)]

Relb_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_wt[,c(856)]
Relb_peaks_wt<-names(Relb_peaks_wt[Relb_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Relb_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Relb_peaks_wt)]

Rela_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(698)]
Rela_peaks_tg<-names(Rela_peaks_tg[Rela_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Rela_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Rela_peaks_tg)]

Rela_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_wt[,c(698)]
Rela_peaks_wt<-names(Rela_peaks_tg[Rela_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Rela_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Rela_peaks_wt)]

Junb_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(127)]
Junb_peaks_tg<-names(Junb_peaks_tg[Junb_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Junb_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Junb_peaks_tg)]

Junb_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(127)]
Junb_peaks_wt<-names(Junb_peaks_wt[Junb_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Junb_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Junb_peaks_wt)]

Bach1_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(119)]
Bach1_peaks_tg<-names(Bach1_peaks_tg[Bach1_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Bach1_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Bach1_peaks_tg)]

Bach1_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_wt[,c(119)]
Bach1_peaks_wt<-names(Bach1_peaks_wt[Bach1_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Bach1_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Bach1_peaks_wt)]

Nfe2l2_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(101)]
Nfe2l2_peaks_tg<-names(Nfe2l2_peaks_tg[Nfe2l2_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Nfe2l2_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Nfe2l2_peaks_tg)]

Nfe2l2_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_wt[,c(101)]
Nfe2l2_peaks_wt<-names(Nfe2l2_peaks_wt[Nfe2l2_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Nfe2l2_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Nfe2l2_peaks_wt)]

Batf3_peaks_tg<-motifMatrix_fibro6a_Co_NEW_mir_tg[,c(789)]
Batf3_peaks_tg<-names(Batf3_peaks_tg[Batf3_peaks_tg==1])
peaks_fibro6a_Co_NEW_mir_tg_Batf3_peaks_tg_all<-peaks_fibro6a_Co_NEW_mir_tg[(elementMetadata(peaks_fibro6a_Co_NEW_mir_tg)[, "peakName"] %in% Batf3_peaks_tg)]

Batf3_peaks_wt<-motifMatrix_fibro6a_Co_NEW_mir_wt[,c(789)]
Batf3_peaks_wt<-names(Batf3_peaks_wt[Batf3_peaks_wt==1])
peaks_fibro6a_Co_NEW_mir_wt_Batf3_peaks_wt_all<-peaks_fibro6a_Co_NEW_mir_wt[(elementMetadata(peaks_fibro6a_Co_NEW_mir_wt)[, "peakName"] %in% Batf3_peaks_wt)]

features_tg_fani_Feb_2022<-GRangesList(Nfe2l2 = peaks_fibro6a_Co_NEW_mir_tg_Nfe2l2_peaks_tg_all, Batf3 = peaks_fibro6a_Co_NEW_mir_tg_Batf3_peaks_tg_all, Junb = peaks_fibro6a_Co_NEW_mir_tg_Junb_peaks_tg_all, Rela = peaks_fibro6a_Co_NEW_mir_tg_Rela_peaks_tg_all, Relb = peaks_fibro6a_Co_NEW_mir_tg_Relb_peaks_tg_all, Bach1 = peaks_fibro6a_Co_NEW_mir_tg_Bach1_peaks_tg_all, Peaks = getPeakSet(fibro6a_Co_NEW_mir_tg))
features_wt_fani_Feb_2022<-GRangesList(Nfe2l2 = peaks_fibro6a_Co_NEW_mir_wt_Nfe2l2_peaks_wt_all, Batf3 = peaks_fibro6a_Co_NEW_mir_wt_Batf3_peaks_wt_all, Junb = peaks_fibro6a_Co_NEW_mir_wt_Junb_peaks_wt_all, Rela = peaks_fibro6a_Co_NEW_mir_wt_Rela_peaks_wt_all, Relb = peaks_fibro6a_Co_NEW_mir_wt_Relb_peaks_wt_all, Bach1 = peaks_fibro6a_Co_NEW_mir_wt_Bach1_peaks_wt_all, Peaks = getPeakSet(fibro6a_Co_NEW_mir_wt))

## Plot Browser Tracks TG ##

pal <- paletteDiscrete(values = fibro6a_Co_NEW_mir_tg$super_clusters)
pooled_colors_new<- c("aquamarine3", "purple4", "firebrick3", "deepskyblue2")
names(pooled_colors_new)<-names(pal)

p221 <- plotBrowserTrack(
  ArchRProj = fibro6a_Co_NEW_mir_tg, 
  groupBy = "super_clusters", 
  geneSymbol = c("Mir221"), upstream = 65000,
  downstream = 65000,
  loops = getPeak2GeneLinks(fibro6a_Co_NEW_mir_tg), pal = pooled_colors_new, features = features_tg_fani_Feb_2022, minCells = 5)
plotPDF(plotList = p221, name = "Mir221-Peak2GeneLinks-tg-Fani-2022-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_tg, 
        addDOC = FALSE, width = 5, height = 3, geneAnnotation=NULL)
# addDOC = FALSE, width = 5, height = 25)
plotPDF(plotList = p221, name = "Mir221-Peak2GeneLinks-tg-Fani-2022-extended-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_tg, 
        addDOC = FALSE, width = 5, height = 12, geneAnnotation=NULL)

p222 <- plotBrowserTrack(
  ArchRProj = fibro6a_Co_NEW_mir_tg, 
  groupBy = "super_clusters", 
  geneSymbol = c("Mir222"), upstream = 65000,
  downstream = 65000,
  loops = getPeak2GeneLinks(fibro6a_Co_NEW_mir_tg), pal = pooled_colors_new, features = features_tg_fani_Feb_2022, minCells = 5)
plotPDF(plotList = p222, name = "Mir222-Peak2GeneLinks-tg-Fani-2022-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_tg, 
        addDOC = FALSE, width = 5, height = 3, geneAnnotation=NULL)
plotPDF(plotList = p222, name = "Mir222-Peak2GeneLinks-tg-Fani-2022-extended-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_tg, 
        addDOC = FALSE, width = 5, height = 12, geneAnnotation=NULL)

## Plot Browser Tracks WT ##

pal <- paletteDiscrete(values = fibro6a_Co_NEW_mir_wt$super_clusters)
pooled_colors_new<- c("aquamarine3", "purple4", "firebrick3", "deepskyblue2")
names(pooled_colors_new)<-names(pal)

p221 <- plotBrowserTrack(
  ArchRProj = fibro6a_Co_NEW_mir_wt, 
  groupBy = "super_clusters", 
  geneSymbol = c("Mir221"), upstream = 65000,
  downstream = 65000,
  loops = getPeak2GeneLinks(fibro6a_Co_NEW_mir_wt), pal = pooled_colors_new, features = features_wt_fani_Feb_2022, minCells = 5)
plotPDF(plotList = p221, name = "Mir221-Peak2GeneLinks-wt-Fani-2022-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_wt, 
        addDOC = FALSE, width = 5, height = 3, geneAnnotation=NULL)
# addDOC = FALSE, width = 5, height = 25)
plotPDF(plotList = p221, name = "Mir221-Peak2GeneLinks-wt-Fani-2022-extended-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_wt, 
        addDOC = FALSE, width = 5, height = 12, geneAnnotation=NULL)

p222 <- plotBrowserTrack(
  ArchRProj = fibro6a_Co_NEW_mir_wt, 
  groupBy = "super_clusters", 
  geneSymbol = c("Mir222"), upstream = 65000,
  downstream = 65000,
  loops = getPeak2GeneLinks(fibro6a_Co_NEW_mir_wt), pal = pooled_colors_new, features = features_wt_fani_Feb_2022, minCells = 5)
plotPDF(plotList = p222, name = "Mir222-Peak2GeneLinks-wt-Fani-2022-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_wt, 
        addDOC = FALSE, width = 5, height = 3, geneAnnotation=NULL)
plotPDF(plotList = p222, name = "Mir222-Peak2GeneLinks-wt-Fani-2022-extended-official.pdf", 
        ArchRProj = fibro6a_Co_NEW_mir_wt, 
        addDOC = FALSE, width = 5, height = 12, geneAnnotation=NULL)

############################


