library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(sctransform)
library(SeuratObject)
library(harmony)
#library(presto)
#options(future.globals.maxSize = 1000 * 1024^10)


Read_10X<-function(path, name){
  data <- Read10X(data.dir = path,gene.column=1)
  data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 2000) 
  data$Time <- name

  data}

Seed4_1<-Read_10X("DAP4_1","Seed4_1")
Seed4_2<-Read_10X("DAP4_2","Seed4_2")
Seed8_1<-Read_10X("DAP8_1","Seed8_1")
Seed8_2<-Read_10X("DAP8_2","Seed8_2")

### Merge
seed_combined_merged <- merge(Seed4_1, y = c( Seed4_2, Seed8_1, Seed8_2))

save(seed_combined_merged, file = "seed.only_merged0331.RData")

seed_combined_merged <- SCTransform(seed_combined_merged)
seed_combined_merged <- RunPCA(object = seed_combined_merged, npcs = 30, verbose = FALSE)
seed_combined_merged <- RunUMAP(object = seed_combined_merged, reduction = "pca", dims = 1:30)
seed_combined_merged <- FindNeighbors(object = seed_combined_merged, reduction = "pca", dims = 1:30)
seed_combined_merged <- FindClusters(seed_combined_merged, resolution = 0.5)  # if change resolution here, make sure you change it down where you transfer identities!
#DimPlot(object = seed_combined_merged, reduction = "umap", label = TRUE,raster=FALSE)

save(seed_combined_merged, file = "seed_combined_merged0331.RData")


###Integration
ifnb.list <- SplitObject(seed_combined_merged,split.by="Time")
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, verbose = FALSE)
}else{
  ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, verbose = FALSE, method = "glmGamPoi")
}

features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
plant_anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:30, normalization.method = "SCT",anchor.features = features)

seed_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30, normalization.method = "SCT")

# # remove(plant_anchors)
# # remove(Seed4_1, Seed4_2, Seed8_1, Seed8_2)
#
DefaultAssay(object = seed_combined_integrated) <- "integrated"
seed_combined_integrated <- ScaleData(object = seed_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
seed_combined_integrated <- RunPCA(object = seed_combined_integrated, npcs = 30, verbose = FALSE)
seed_combined_integrated <- RunUMAP(object = seed_combined_integrated, reduction = "pca", dims = 1:30)
seed_combined_integrated <- FindNeighbors(object = seed_combined_integrated, reduction = "pca", dims = 1:30)
seed_combined_integrated <- FindClusters(seed_combined_integrated, resolution = 0.5)  # if change resolution here, make sure you change it down where you transfer identities!
save(seed_combined_integrated, file = "seed_combined_integrated0331.RData")
pdf("seed.cluster.pdf")
DimPlot(object = seed_combined_integrated, reduction = "umap", label = TRUE,raster=FALSE)
dev.off()
load("seed_combined_integrated0331.RData")

#seed_combined_integrated <-PrepSCTFindMarkers(seed_combined_integrated)
DefaultAssay(seed_combined_integrated)
seed_combined_integrated <-PrepSCTFindMarkers(object = seed_combined_integrated)
marker=FindAllMarkers(seed_combined_integrated,assay="SCT",only.pos=T,min.pct = 0.1,logfc.threshold=0.25)
write_tsv(marker,"integration.SCT.filter.Marker.resolution300.5.txt")

marker=FindAllMarkers(seed_combined_integrated,assay="RNA",only.pos=T,min.pct = 0.1,logfc.threshold=0.25)
write_tsv(marker,"integration.RNA.filter.Marker.resolution300.5.txt")
