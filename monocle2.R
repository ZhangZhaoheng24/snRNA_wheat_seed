rm(list = ls())
library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)

load("seed_combined_integrated0331.RData")
table(Idents(seed_combined_integrated))
seed_combined_integrated@meta.data
table(seed_combined_integrated$orig.ident)
head(seed_combined_integrated@meta.data)
#scRNA=seed_combined_integrated
#rm(seed_combined_integrated)
#
levels(Idents(seed_combined_integrated)) 
all_cells <- colnames(seed_combined_integrated)
set.seed(123)  
selected_cells <- sample(all_cells, size = floor(length(all_cells) * 0.5), replace = FALSE)
scRNA <- subset(seed_combined_integrated, cells = selected_cells)
scRNA$celltype = Idents(scRNA)
table(Idents(scRNA))

## 
## CD8+ T-cells     NK cells 
##          544           92
set.seed(1234)

ct <- scRNA@assays$RNA@counts

gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)

pd <- new("AnnotatedDataFrame",
          data=scRNA@meta.data)

sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds


sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)


dd = "differentialGeneTestall0416.Rdata"


diff_test_res <- differentialGeneTest(sc_cds,
                                     fullModelFormulaStr = " ~ celltype + Time",
                                     reducedModelFormulaStr = " ~ Time",
                                     cores = 32)
save(diff_test_res,file = dd)

load(dd)


ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))


head(ordering_genes)
length(ordering_genes)


sc_cds <- setOrderingFilter(sc_cds, ordering_genes)

plot_ordering_genes(sc_cds)
sc_cds <- reduceDimension(sc_cds,residualModelFormulaStr = "~Time")
sc_cds <- orderCells(sc_cds)
save(sc_cds,file = "cluster_all_psutime0416.RData")