# Continued from 010_QC.R
analysis_step <- "020_clustering"

# load packages ----
library(tidyverse)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
seu_all <- readRDS(file.path("RDSfiles", "seu_010_filt.RDS"))
load(file.path("RDSfiles", "cellgroup_names.Rdata"))

# Epithelial ----
# Make directories
plot_path <- file.path("plot", analysis_step, "epi")
res_path <- file.path("result", analysis_step, "epi")
fs::dir_create(c(plot_path, res_path))

# Clustering
seu <- seu_all[, epi_names]
seu$cellgroup <- "Epi."
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)    # adjust resolution when necessary
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3.5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "colon_epi_markers.txt"))

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
         width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
markers <- FindAllMarkers(seu, only.pos = TRUE)
# add_feat <- "Fabp6"
# FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# ggsave(paste0(add_feat, ".png"), path = fp_path, width = 5, height = 5, units = "in", dpi = 150)
# rm(add_feat)
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Add celltype annotation and save the Seurat object
seu <- subset(seu, idents = 15, invert = TRUE)
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(0,2,3,12)] <- "Tum.1"
seu$celltype[seu$seurat_clusters %in% c(1)] <- "Tum.2"
seu$celltype[seu$seurat_clusters %in% c(4,7,10)] <- "Tum.3"
seu$celltype[seu$seurat_clusters %in% c(5)] <- "Tum.4"
seu$celltype[seu$seurat_clusters %in% c(8)] <- "TA"    # Dmbt1, Aqp4, partly Mki67
seu$celltype[seu$seurat_clusters %in% c(11)] <- "Colono"
seu$celltype[seu$seurat_clusters %in% c(14)] <- "Sec.Pro."    # Muc2, Mki67
seu$celltype[seu$seurat_clusters %in% c(6,13)] <- "Goblet"
seu$celltype[seu$seurat_clusters %in% c(9)] <- "Paneth"
seu$celltype[seu$seurat_clusters %in% c(16)] <- "EEC"
# seu$celltype[seu$seurat_clusters %in% c(15)] <- "Pancreas"
seu$celltype <- factor(seu$celltype, levels = c("Tum.1", "Tum.2", "Tum.3", "Tum.4", 
                                                "TA", "Colono", "Sec.Pro.", "Goblet",
                                                "Paneth", "EEC"))
DimPlot(seu, group.by = "celltype", cols = "alphabet", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)

# Dot plot
features <- readLines(file.path("aux_data", "gene_set", "colon_epi_markers_dp.txt"))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))

# Stromal ----
# Make directories
plot_path <- file.path("plot", analysis_step, "str")
res_path <- file.path("result", analysis_step, "str")
fs::dir_create(c(plot_path, res_path))

# Clustering
seu <- seu_all[, str_names]
seu$cellgroup <- "Str."
npcs <- 15    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)    # adjust resolution when necessary
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "colon_str_markers.txt"))

# save_fp <- function(feature, seu, path){
#   tryCatch({
#     p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
#       NoAxes() + NoLegend()
#     ggsave(paste0(feature, ".png"), plot = p, path = path, 
#            width = 5, height = 5, units = "in", dpi = 150)
#   }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
# }
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
markers <- FindAllMarkers(seu, only.pos = TRUE)
# add_feat <- "Myh11"
# FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# ggsave(paste0(add_feat, ".png"), path = fp_path, width = 5, height = 5, units = "in", dpi = 150)
# rm(add_feat)
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Add celltype annotation and save the Seurat object
seu <- subset(seu, idents = c(0,11,14), invert = TRUE)    # These clusters seem epithelial or immune
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(7)] <- "Fib.1"
seu$celltype[seu$seurat_clusters %in% c(1,2,4)] <- "Fib.2"
seu$celltype[seu$seurat_clusters %in% c(8)] <- "Fib.3"
seu$celltype[seu$seurat_clusters %in% c(12)] <- "Prolif."
seu$celltype[seu$seurat_clusters %in% c(3,9)] <- "Myo."
seu$celltype[seu$seurat_clusters %in% c(6)] <- "BEC"
seu$celltype[seu$seurat_clusters %in% c(5)] <- "LEC"
seu$celltype[seu$seurat_clusters %in% c(10)] <- "ICC"
seu$celltype[seu$seurat_clusters %in% c(13)] <- "Meso."
seu$celltype <- factor(seu$celltype, levels = c("Fib.1", "Fib.2", "Fib.3", "Prolif.", 
                                                "Myo.", "BEC", "LEC", "ICC",
                                                "Meso."))
DimPlot(seu, group.by = "celltype", cols = "alphabet", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)

# Dot plot
features <- readLines(file.path("aux_data", "gene_set", "colon_str_markers_dp.txt"))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_020_str.RDS"))

# Immune ----
# Make directories
plot_path <- file.path("plot", analysis_step, "imm")
res_path <- file.path("result", analysis_step, "imm")
fs::dir_create(c(plot_path, res_path))

# Clustering
seu <- seu_all[, imm_names]
seu$cellgroup <- "Imm."
npcs <- 20    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)    # adjust resolution when necessary
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "alphabet") + NoAxes()
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("sample.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines(file.path("aux_data", "gene_set", "colon_imm_markers.txt"))

# save_fp <- function(feature, seu, path){
#   tryCatch({
#     p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
#       NoAxes() + NoLegend()
#     ggsave(paste0(feature, ".png"), plot = p, path = path, 
#            width = 5, height = 5, units = "in", dpi = 150)
#   }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
# }
sapply(features, save_fp, seu, fp_path)

# Check markers interactively when necessary
markers <- FindAllMarkers(seu, only.pos = TRUE)
add_feat <- "Ms4a2"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 5, height = 5, units = "in", dpi = 150)
rm(add_feat)
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Add celltype annotation and save the Seurat object
seu <- subset(seu, idents = c(5,10), invert = TRUE)    # These clusters seem epithelial
seu$celltype <- ""
seu$celltype[seu$seurat_clusters %in% c(0,7,8)] <- "Bcell"
seu$celltype[seu$seurat_clusters %in% c(4)] <- "Plasma"
seu$celltype[seu$seurat_clusters %in% c(2,3)] <- "CD4-T"
seu$celltype[seu$seurat_clusters %in% c(9)] <- "CD8-T"
seu$celltype[seu$seurat_clusters %in% c(16)] <- "Prolif."
seu$celltype[seu$seurat_clusters %in% c(1,12,13)] <- "Macro."
seu$celltype[seu$seurat_clusters %in% c(11)] <- "Neutro."
seu$celltype[seu$seurat_clusters %in% c(14,15)] <- "DC"
seu$celltype[seu$seurat_clusters %in% c(6)] <- "Mast"
seu$celltype <- factor(seu$celltype, levels = c("Bcell", "Plasma", "CD4-T", "CD8-T", 
                                                "Prolif.", "Macro.", "Neutro.", "DC",
                                                "Mast"))
DimPlot(seu, group.by = "celltype", cols = "alphabet", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)

# Dot plot
features <- readLines(file.path("aux_data", "gene_set", "colon_imm_markers_dp.txt"))
DotPlot(seu, group.by = "celltype", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 6.5, height = 4, units = "in", dpi = 150)

saveRDS(seu, file = file.path("RDSfiles", "seu_020_imm.RDS"))
