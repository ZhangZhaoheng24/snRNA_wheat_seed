library(Seurat)
library(tidyverse)
library(cowplot)
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(harmony)
set.seed(12345)
setwd("/home/zhangzh/en_sc")
remove(scRNA_harmony)
load("seed_combined_integrated_UcellbulkRNA_ATAC_jisu.RData")
load(file="seed_hdWGCNA0604seed_4_n.RData")
remove(sce,scRNA,scRNA_harmony,all)
DefaultAssay(object = seed_combined_integrated) <-"RNA"

seurat_obj <- SetupForWGCNA(seed_combined_integrated,
                            gene_select = "fraction",
                            fraction = 0.05,
                            wgcna_name = "Metastates")
wgcna_name="Metastates"
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_types", "Time2"),
  reduction = 'umap',
  k = 25,
  max_shared = 10,
  ident.group = 'cell_types'
)
table(seurat_obj$Time2, seurat_obj$metacell_grouping)

seurat_obj <- NormalizeMetacells(seurat_obj)

rm(seed_combined_integrated)

rm(metacell_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  #group_name = "BM",
  #group.by = 'group',
  assay = 'RNA',
  layer = 'data'
)
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)


library(patchwork)
# plot_list <- PlotSoftPowers(seurat_obj)
# wrap_plots(plot_list, ncol = 2)
# ggsave("hdWGCNA_soft_power0604h.pdf")

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 4,  
  setDatExpr = FALSE,
  tom_name = 'seed_4_n',
  networkType = "unsigned"
)
modules <- GetModules(seurat_obj, "Metastates")
net <- GetNetworkData(seurat_obj, "Metastates")
head(modules)
head(seurat_obj)
order=net$dendrograms[[1]]$order
modules2=modules[order,]
write.table(modules2,"/home/zhangzh/en_sc/modules_seed_4_n.txt",sep='\t',quote=F)
head(modules2)
aveexp=read.csv("all_AggregateExpression_c_t.txt",header=T)
rownames(aveexp)=aveexp$X
aveexp=aveexp[,-1]
head(aveexp)
#aveexp0=log2(aveexp+1)
aveexp0=t(scale(t(aveexp), center = TRUE, scale = TRUE))
#head(aveexp0)
#aveexp2=t(aveexp[modules2$gene_name,])
aveexp2=t(aveexp0[modules2$gene_name,])
aveexp3=aveexp2
#aveexp3=t(scale(t(aveexp2), center = TRUE, scale = TRUE))
head(aveexp3[1:10,1:10])
cols=as.data.frame(rownames(aveexp2))
colnames(cols)="ID"
library(tidyr)
cols2 <- cols %>%
  separate(
    ID,              
    into = c("celltype", "Time"),  
    sep = "_",              
    remove = FALSE        
  )

head(aveexp2[1:10,1:10])
library(ComplexHeatmap)
rownames(cols2)=cols2$ID
cols2=cols2[order(cols2$Time),]
cols2=cols2[order(cols2$celltype),]
aveexp3=aveexp3[cols2$ID,]
head(aveexp3[1:10,1:10])
col_fun = circlize::colorRamp2(c(0, 2), c("white", "#FF0504"))
p1 <- Heatmap(aveexp3, name = "Z-score", 
              show_row_names = F, show_column_names = F,
              #row_split=cols2[rownames(aveexp3), c("celltype","Time")],
              right_annotation = rowAnnotation(df = cols2[,c("celltype","Time"),drop=F],
              col = list(
              celltype=c("RNA.AL"="#8A7765","RNA.Central"="#68AEE0","RNA.CF"="#B8D168","RNA.Dividing.cell"="#188645",
                         "RNA.Embryo"="#EA7368","RNA.ESR"="#DB7FB1","RNA.Pericarp"="#206264","RNA.SC"="#694E83",
                         "RNA.Scutellum"="#DF9339","RNA.TCSE"="#F4CF5A","RNA.Unknown"="#E5E5E4"),
              Time=c("Seed4"="#599799","Seed8"="#E99573")
              ),simple_anno_size = unit(2, "mm")
              ),
              #column_split = column_groups,
              # cluster_row_slices = TRUE,      # 每个slice内做聚类
              # clustering_method_rows = "complete",
              use_raster = TRUE,
              
              col = col_fun, 
              cluster_columns = FALSE, cluster_rows = FALSE,
              border = FALSE, # border_gp = gpar(col = "black", lwd = 2), #设置热图的边框与粗度
              
              column_names_rot = 45,
              #right_annotation = row_anno,
              column_title_gp =gpar(fontsize = 3),
              rect_gp = gpar(col = NA),# 中间线
              column_names_gp = gpar(fontsize = 14),
              width = unit(15,'cm'),
              height =  unit(4,'cm'),gap = unit(0.5,'mm'))
pdf("heatmap_sc_dram.pdf",width=15,height=5)
p1
dev.off()


PlotDendrogram(seurat_obj, main = 'Seed hdWGCNA Dendrogram')
seurat_obj@misc$active_wgcna
ggsave("hdWGCNA_Dendrogram0604hseed_4_n.pdf")

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "Time2"
)


seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_types', 
  group_name = as.vector(unique(seurat_obj$cell_types))
)


PlotKMEs(seurat_obj, ncol = 5, n_hubs = 10)


modules <- GetModules(seurat_obj)


hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 50)
write.table(hub_df,"hub_gene_in_net50seed_4_n.txt",sep='\t',quote=F)
modules=modules[modules$module!="grey",]
write.table(modules,"modules_gene_in_netnograyseed_4_n.txt",sep='\t',quote=F)
#
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 50,  
  method = "Seurat"  
)


plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs',
  order = TRUE
)
wrap_plots(plot_list, ncol=5)
ggsave("hdWGCNA_hubgene_scoreseed_4_n.pdf",width = 15,height = 7)

seurat_obj@active.ident <- seurat_obj$cell_types


wrap_plots(plot_list, ncol = 6)
ModuleCorrelogram(seurat_obj)


MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']


seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# 绘制 DotPlot
p <- DotPlot(
  seurat_obj,
  features = mods,
  group.by = "cell_types"
) +
  RotatedAxis() +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")

p
ggsave("hdWGCNA_hubgene_score_dotplotseed_4_n.pdf",width = 9,height = 4)

# p <- VlnPlot(
#   seurat_obj,
#   features = "black",
#   group.by = "cell_types",
#   pt.size = 0
# ) +
#   geom_boxplot(width = .25, fill = "white") +
#   xlab("") +
#   ylab("hME") +
#   NoLegend()
# 
# p
# save(seurat_obj,file="seed_hdWGCNA0604seed_4_n.RData")
for (mod in unique(hub_df$module)){
  ModuleNetworkPlot(seurat_obj = seurat_obj,   
                    n_inner = 15,
                    n_outer = 17,
                    n_conns = 500,
                    plot_size = c(6, 6),
                    mods = mod,outdir = "ModuleNetworksseed_4_n")
}

load("seed_hdWGCNA0604seed_4_n.RData")

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
umap_df <- GetModuleUMAP(seurat_obj)
write.table(umap_df,"seed_hdWGCNA0604seed_4_n_umap.txt",sep='\t',quote = F)
# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
ggsave("hdWGCNA_UMAP_networkSeed_4_n10.pdf",width = 8,height = 7)
options(future.globals.maxSize = 50442450944) 

head(umap_df)
pdf("hdWGCNA_UMAP_network2_Seed_4_n10_1218.pdf",width = 6,height = 15)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.05, 
  label_hubs=0 ,#
  keep_grey_edges=FALSE,
  label_genes=c("TraesCS1A02G340400","TraesCS1B02G273300","TraesCS2B02G194200","TraesCS2A02G385900","TraesCS5A02G422900","TraesCS2B02G530500","TraesCS5A02G472000","TraesCS3B02G363200","TraesCS1A02G038200","TraesCS3A02G020300","TraesCS3A02G344700")
  
)


dev.off()


pdf("hdWGCNA_hub_network2Seed_4_n.pdf",width = 8,height = 7)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 5,       # 
  n_other = 15,     # 
  edge_prop = 0.75, # 
  mods = "all"      #
)
dev.off()
