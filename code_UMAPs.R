# Generating Seurat objects and UMAPs 

library(Seurat)
library(ggplot2)

# all cells =================
# load csv files
count_matrix <- read.csv("count_matrix.csv", row.names = 1) 
metadata <- read.csv("metadata.csv", row.names = 1)  

# assemble Seurat object
seurat_obj <- CreateSeuratObject(counts = count_matrix, meta.data = metadata)

# plot publication UMAP
umap_coords <- read.csv("umap_embeddings_all_cells.csv", row.names = 1)  # Ensure row names are cell names to match the Seurat object
rownames(umap_coords) <- gsub("-", ".", rownames(umap_coords))
seurat_obj@reductions$umap <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords), assay = "RNA", key = "UMAP_")

plot <- DimPlot(seurat_obj, group.by = "annotation_publication", reduction = "umap",
                cols = c("#f9a81d", # Astrocytes 
                         "#feaf8a", # BergmannGlia 
                         "#fd7a8c", # ChoroidPlexus
                         "#cc89d6", # Ependymal
                         "#bfcff0", # Leptomeningeal 
                         "#9ce7c9", # Macrophages
                         "#4dc656", # Microglia 
                         "#a6aab2", # Oligo
                         "#6098b4", # Peripheral Immune cells
                         "#dbb7af" # Vascular
                ), raster = F
) + theme_void()

plot

# microglia ============

# subset Seurat object
Idents(seurat_obj) = "annotation_publication"
microglia = subset(seurat_obj, idents = "Microglia")

# plot publication UMAP
umap_coords <- read.csv("umap_embeddings_microglia.csv", row.names = 1) 
rownames(umap_coords) <- gsub("-", ".", rownames(umap_coords))
microglia@reductions$umap <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords), assay = "RNA", key = "UMAP_")


plot <- DimPlot(microglia, group.by = "annotation_publication_microglia", reduction = "umap",
                cols = c("#9bcfe4", # cluster1 
                         "#a6e07f", # cluster2
                         "#007ab9", # cluster3
                         "#f70000", # cluster4
                         "#ffbb5e", # cluster5
                         "#f49697", # cluster6
                         "#00a200" # cluster7
                ), raster = F
) + theme_void()

plot
