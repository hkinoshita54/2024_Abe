analysis_step <- "010_QC"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)
library(DoubletFinder)

# Make directories ----
fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
## sample_name should be identical to the folder name of each sample
sample_table <- read_excel(file.path("data", "sample_table.xlsx")) %>% as.data.frame
sample_name <- sample_table$sample_name
nsample <- length(sample_name)
exp_group <- sample_table$exp_group
rm(sample_table)

ReadData <- function(sample_name){
  data_dir <- file.path("data", sample_name)
  cts <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(counts = cts, project = sample_name, min.cells = 3, min.features = 200)
  return(seu)
}
seu_list <-lapply(sample_name, ReadData)

# DoubletFinder ----
source(file.path("Rscripts", "011_DoubletFinder.R"))

## merge the list
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)], add.cell.ids = sample_name)
rm(seu_list)
seu@meta.data <- seu@meta.data[,c(1:3,7)]

# QC & filter ----
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_unfil.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

seu <- subset(seu, subset = doublet_finder == "Singlet")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_singlet.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_filt.png", path = plot_path, width = 3*nsample, height = 3, units = "in", dpi = 150)

# Save filtered object ----
seu <- JoinLayers(seu)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
saveRDS(seu, file.path("RDSfiles", "seu_010_filt.RDS"))
seu_all<-seu

# Initial Clustering ----
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots ----
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "colon_group_markers.txt"))

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
         width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary ----
# markers <- FindAllMarkers(seu, only.pos = TRUE)
# add_feat <- "Rbpjl"
# FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# ggsave(paste0(add_feat, ".png"), path = fp_path, width = 5, height = 5, units = "in", dpi = 150)
# rm(add_feat)

# Add cellgroup annotation ----
epi_names <- colnames(seu)[seu$seurat_clusters %in% c(0:2,4,5,9:12,16,18)]
str_names <- colnames(seu)[seu$seurat_clusters %in% c(3,13,19,20)]
imm_names <- colnames(seu)[seu$seurat_clusters %in% c(6:8,14,15,17)]
save(epi_names, str_names, imm_names, file = file.path("RDSfiles", "cellgroup_names.Rdata"))

seu$cellgroup <- "Epi."
seu$cellgroup[str_names] <- "Str."
seu$cellgroup[imm_names] <- "Imm."
seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Str.", "Imm."))
DimPlot(seu, group.by = "cellgroup") & NoAxes()
ggsave("cellgroup.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
features <- readLines(file.path("aux_data", "gene_set", "colon_group_markers_dp.txt"))
DotPlot(seu, group.by = "cellgroup", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)
