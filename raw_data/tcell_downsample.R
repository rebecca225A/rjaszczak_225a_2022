library(Seurat)
library(dittoSeq)
library(dplyr)

# Users running this code must update this directory 
download_dir = "/Users/rebecca/Downloads/"

# In this directory you'll need:
# 1)
# Dan's annotated Tcells.rds
# https://figshare.com/articles/dataset/Processed_naive_T_cell_single-cell_RNA-seq_Seurat_object/11886891
# 2)
# Raw cellranger outs and demux info UNZIPPED:
# https://figshare.com/articles/dataset/naive_T_cell_single-cell_RNA-seq_raw_counts_and_annotation/11894637
# Tcells should be unzipped DIRECTLY into download_dir.
# you'll have to modify the paths if your system unzips in a nested manner, eg: Tcells/Tcells/

# read in Dan's analyzed rds
Tcells_analyzed = readRDS(file.path(download_dir,"Tcells.rds"))
Idents(Tcells_analyzed) = "Sample"

# quantify number of cells in each identiy and create downsampling values
tcell_table = as.data.frame(table(Tcells_analyzed@meta.data$Sample)) %>% 
  dplyr::rename(Freq_Full = Freq) %>%
  dplyr::mutate(Freq_10 = round(Freq_Full/10)) %>%
  dplyr::rename(Sample = Var1)

# make subset cell lists
set.seed(2598) #set seed for re-producability 

Tcells_analyzed #40152 cells
40152/10 #4015.2 target

# create an empty list and add downsampled cell IDs to it from each identity according to the amounts previously determined
tcell_subset_10_list = c()
for ( x in 1:length(tcell_table$Sample) ) {
  cells_ident = WhichCells(Tcells_analyzed, idents = tcell_table$Sample[x], downsample = tcell_table$Freq_10[x])
  tcell_subset_10_list = c(tcell_subset_10_list, cells_ident)
}

#subset analyzed object
Tcells_analyzed_ds10 = Tcells_analyzed[,tcell_subset_10_list]

#compare full vs downsampled object
#this might take a hot second to generate depending on how much RAM your laptop has
#full = dittoDimPlot(Tcells_analyzed, "Sample", split.by = "Sample", split.ncol = 5)
#ds = dittoDimPlot(Tcells_analyzed_ds10, "Sample", split.by = "Sample", split.ncol = 5)
#full + ds 

# initialize a raw object with "bad" quality cells, for qc practice
Tcells_raw <- CreateSeuratObject(Read10X(file.path(download_dir,"Tcells/cellranger_Raw/")))

# demux raw object
Tcells_demuxed <- importDemux(
  Tcells_raw,
  demuxlet.best = c(file.path(download_dir,"Tcells/Demuxlet/CD4.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD4-8.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD8.best")),
  lane.names = c("CD4","CD4-8","CD8"))

#Remove the "CD4_" at the start of sample names (due to coming from RNAseq data with those names)
Tcells_demuxed[["Sample"]] <- sapply(
  meta("Sample",Tcells_demuxed),
  function(X) strsplit(X, split = "CD4_")[[1]][2])

#summarize snp data
demux.SNP.summary(
  Tcells_demuxed,
  plots = c("jitter","vlnplot","boxplot"),
  boxplot.color = "white",
  boxplot.fill = F,
  add.line = 50)
demux.calls.summary(Tcells_demuxed, singlets.only = FALSE)

#subset based on demuxlet singlets
Tcells_demuxed <- subset(Tcells_demuxed, subset = demux.doublet.call == "SNG")
#this object contains 42793 cells (singlets, but may be poorer quality)

# make a highly qc'd object
Tcells_raw[["percent.mito"]] <- PercentageFeatureSet(Tcells_raw, pattern = "^MT-")
Tcells_raw[["percent.ribo"]] <- PercentageFeatureSet(Tcells_raw, pattern = "^RPS|^RPL")

Tcells.cut <- subset(Tcells_raw, subset = nFeature_RNA > 750)
Tcells.cut <- subset(Tcells.cut, subset = nCount_RNA > 1500)
Tcells.cut <- subset(Tcells.cut, subset = percent.mito < 5)
Tcells.cut
Tcells.cut <- importDemux(
  Tcells.cut,
  demuxlet.best = c(file.path(download_dir,"Tcells/Demuxlet/CD4.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD4-8.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD8.best")),
  lane.names = c("CD4","CD4-8","CD8"))

Tcells.cut[["Sample"]] <- sapply(
  meta("Sample",Tcells.cut),
  function(X) strsplit(X, split = "CD4_")[[1]][2])

demux.SNP.summary(
  Tcells.cut,
  plots = c("jitter","vlnplot","boxplot"),
  boxplot.color = "white",
  boxplot.fill = F,
  add.line = 50)
demux.calls.summary(Tcells.cut, singlets.only = FALSE)
Tcells.cut <- subset(Tcells.cut, subset = demux.doublet.call == "SNG")

#create vectors of demuxed vs "cut" or qc'd cells
demuxed_cells_no_qc = Cells(Tcells_demuxed)
demuxed_cells_cut = Cells(Tcells.cut)

# sanity check that all cut cells are in the unfiltered, no qc object
table(demuxed_cells_cut %in% demuxed_cells_no_qc)

# cell ids that will fail the filtering based on nFeature_RNA, etc
fail_qc = demuxed_cells_no_qc[!demuxed_cells_no_qc %in% demuxed_cells_cut]

# subset 15% of these poor quality cells for inclusion in class object, for filtering practice
poor_quality_cells = sample(demuxed_cells_no_qc, length(fail_qc)*.15)
high_mito = Cells(subset(Tcells_raw, subset = percent.mito > 5)) #tons of theses aren't annotated in demux so taking all

# combine 10% downsample of quality cells and low quality cells
full_subset_list = c(tcell_subset_10_list, poor_quality_cells, high_mito)

# remove duplicates
full_subset_list = full_subset_list[!duplicated(full_subset_list)]
table(duplicated(full_subset_list))

#make a seurat object subset based on 10% downsample of quality cells and low quality cells
Tcells_raw_ds10 = Tcells_raw[,full_subset_list]

#demux ds object
Tcells_raw_ds10 <- importDemux(
  Tcells_raw_ds10,
  demuxlet.best = c(file.path(download_dir,"Tcells/Demuxlet/CD4.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD4-8.best"),
                    file.path(download_dir,"Tcells/Demuxlet/CD8.best")),
  lane.names = c("CD4","CD4-8","CD8"))

#Remove the "CD4_" at the start of sample names (due to coming from RNAseq data with those names)
Tcells_raw_ds10[["Sample"]] <- sapply(
  meta("Sample",Tcells_raw_ds10),
  function(X) strsplit(X, split = "CD4_")[[1]][2])

Tcells_raw_ds10 <- subset(Tcells_raw_ds10, subset = demux.doublet.call == "SNG")

#write out downsampled RDS for class
saveRDS(Tcells_raw_ds10, file = file.path(download_dir,"downsampled_Tcells_raw.rds"))