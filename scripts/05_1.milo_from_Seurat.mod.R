library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)


### 4.3 From Seurat object
library(Seurat)
## Load in the data
alldata <- readRDS(file = "../harmony_seurat/Harmony_round1.rds")

set.seed(42)


## add metadata
# age
alldata$age <- NA
alldata$age[alldata$Sample == "SLE07"] <- 55
alldata$age[alldata$Sample == "SLE01"] <- 57
alldata$age[alldata$Sample == "SLE03"] <- 36
alldata$age[alldata$Sample == "SLE14"] <- 42
alldata$age[alldata$Sample == "SLE11"] <- 37
alldata$age[alldata$Sample == "SLE12"] <- 37
alldata$age[alldata$Sample == "SLE04"] <- 62
alldata$age[alldata$Sample == "SLE13"] <- 80
alldata$age[alldata$Sample == "SLE08"] <- 50
alldata$age[alldata$Sample == "SLE15"] <- 33
alldata$age[alldata$Sample == "SLE06"] <- 61
alldata$age[alldata$Sample == "SLE05"] <- 67
alldata$age[alldata$Sample == "SLE10"] <- 36
alldata$age[alldata$Sample == "SLE02"] <- 42
alldata$age[alldata$Sample == "SLE09"] <- 46
# sex
alldata$sex <- NA
alldata$sex[alldata$Sample == "SLE07"] <- "F"
alldata$sex[alldata$Sample == "SLE01"] <- "F"
alldata$sex[alldata$Sample == "SLE03"] <- "F"
alldata$sex[alldata$Sample == "SLE14"] <- "M"
alldata$sex[alldata$Sample == "SLE11"] <- "F"
alldata$sex[alldata$Sample == "SLE12"] <- "F"
alldata$sex[alldata$Sample == "SLE04"] <- "M"
alldata$sex[alldata$Sample == "SLE13"] <- "M"
alldata$sex[alldata$Sample == "SLE08"] <- "F"
alldata$sex[alldata$Sample == "SLE15"] <- "M"
alldata$sex[alldata$Sample == "SLE06"] <- "F"
alldata$sex[alldata$Sample == "SLE05"] <- "F"
alldata$sex[alldata$Sample == "SLE10"] <- "F"
alldata$sex[alldata$Sample == "SLE02"] <- "F"
alldata$sex[alldata$Sample == "SLE09"] <- "F"

library(dplyr)
alldata$sample_id <- factor(alldata$Sample)
alldata$group_id <- alldata$SampleGroup %>% stringr::str_replace_all(pattern="-", replacement="_") %>% factor()
alldata$institute_id <- factor(alldata$Institute)

alldata$scaled_age <- scale(alldata$age)
alldata
##  An object of class Seurat
##  76712 features across 66915 samples within 7 assays
##  Active assay: RNA (27204 features, 0 variable features)
##   6 other assays present: SCT, refAssay, prediction.score.celltype.l2, prediction.score.celltype.l1, prediction.score.celltype.l3, impADT
##   3 dimensional reductions calculated: pca, harmony, umap

png("SCT.pca.sex.png", width = 800, height = 500)
DimPlot(alldata, reduction = "pca", group.by="sex")
dev.off()

png("SCT.pca.scaled_age.png", width = 800, height = 500)
DimPlot(alldata, reduction = "pca", group.by="scaled_age")
dev.off()

png("Harmony.sex.png", width = 800, height = 500)
DimPlot(alldata, reduction = "harmony", group.by="sex")
dev.off()
png("Harmony.scaled_age.png", width = 800, height = 500)
DimPlot(alldata, reduction = "harmony", group.by="scaled_age")
dev.off()
png("Harmony.institute_id.png", width = 800, height = 500)
DimPlot(alldata, reduction = "harmony", group.by="institute_id")
dev.off()


### 各サンプルか否か
df <- ifelse(alldata$sample_id == "SLE01", "SLE01", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE01")
df <- ifelse(alldata$sample_id == "SLE02", "SLE02", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE02")
df <- ifelse(alldata$sample_id == "SLE03", "SLE03", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE03")
df <- ifelse(alldata$sample_id == "SLE04", "SLE04", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE04")
df <- ifelse(alldata$sample_id == "SLE05", "SLE05", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE05")
df <- ifelse(alldata$sample_id == "SLE06", "SLE06", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE06")
df <- ifelse(alldata$sample_id == "SLE07", "SLE07", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE07")
df <- ifelse(alldata$sample_id == "SLE08", "SLE08", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE08")
df <- ifelse(alldata$sample_id == "SLE09", "SLE09", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE09")
df <- ifelse(alldata$sample_id == "SLE10", "SLE10", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE10")
df <- ifelse(alldata$sample_id == "SLE11", "SLE11", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE11")
df <- ifelse(alldata$sample_id == "SLE12", "SLE12", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE12")
df <- ifelse(alldata$sample_id == "SLE13", "SLE13", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE13")
df <- ifelse(alldata$sample_id == "SLE14", "SLE14", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE14")
df <- ifelse(alldata$sample_id == "SLE15", "SLE15", "OTHER")
alldata <- AddMetaData(object=alldata, metadata=df, col.name="SLE15")





alldata_sce <- as.SingleCellExperiment(alldata)

reducedDim(alldata_sce, "PCA", withDimnames=TRUE) <- alldata[['pca']]@cell.embeddings
reducedDim(alldata_sce, "UMAP", withDimnames=TRUE) <- alldata[['umap']]@cell.embeddings
reducedDim(alldata_sce, "HARMONY", withDimnames=TRUE) <- alldata[['harmony']]@cell.embeddings

set.seed(42)
alldata_milo <- Milo(alldata_sce)

### 5   Construct KNN graph
alldata_milo <- buildGraph(alldata_milo, k = 30, d = 30, reduced.dim = "HARMONY")
##  Constructing kNN graph with k:30


### 6   1. Defining representative neighbourhoods
alldata_milo <- makeNhoods(alldata_milo, prop = 0.05, k = 30, d=30, refined = TRUE, reduced_dims = "HARMONY")
##  Checking valid object

p <- plotNhoodSizeHist(alldata_milo)
ggsave(file="mod/plotNhoodSizeHist.byHARMONY.pdf", plot=p, width=7, height=7)


### 7   Counting cells in neighbourhoods
alldata_milo <- countCells(alldata_milo, meta.data = as.data.frame(colData(alldata_milo)), sample="sample_id")
##  Checking meta.data validity
##  Counting cells in neighbourhoods



### 8   Differential abundance testing
alldata_design <- data.frame(colData(alldata_milo))[,c("sample_id", "scaled_age", "sex", "institute_id", "group_id")]

## Convert batch info from integer to factor
alldata_design$sex <- factor(alldata_design$sex, levels=c("F", "M")) 

alldata_design <- distinct(alldata_design)
rownames(alldata_design) <- alldata_design$sample_id

alldata_design$group_id <- relevel(alldata_design$group_id, "Ctrl")

## Reorder rownames to match columns of nhoodCounts(milo)
alldata_design <- alldata_design[colnames(nhoodCounts(alldata_milo)), , drop=FALSE]

### Computing neighbourhood connectivity
alldata_milo <- calcNhoodDistance(alldata_milo, d=30, reduced.dim = "HARMONY")
##  as(<dgTMatrix>, "dgCMatrix") is deprecated since Matrix 1.5-0; do as(., "CsparseMatrix") instead

### Testing
da_results <- testNhoods(alldata_milo, design = ~ sex + scaled_age + institute_id + group_id, design.df = alldata_design)
##  Using TMM normalisation
##  Performing spatial FDR correction withk-distance weighting

da_results %>%
  arrange(SpatialFDR) %>%
  head() 
##        logFC   logCPM          F       PValue        FDR Nhood SpatialFDR
##  1 6.0200782 8.396771 23.2826525 2.947102e-05 0.08304935   581 0.07758075
##  2 5.8526779 8.203143 18.3242737 2.287877e-04 0.21912442  2423 0.23239556
##  3 5.1873512 8.351078 16.9600717 2.332765e-04 0.21912442  2437 0.23239556
##  4 4.2740054 8.490761 14.5313319 5.585653e-04 0.39350926  2417 0.41199503
##  5 0.9159934 8.241663  0.7474224 3.934060e-01 0.99970042     1 0.99970042
##  6 0.6391230 8.661902  0.2407085 6.268610e-01 0.99970042     2 0.99970042



### Inspecting DA testing results
p <- ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggsave(file="mod/PValue_histo.pdf", plot=p, width=7, height=7)

da2_results <- testNhoods(alldata_milo, design = ~ sex + scaled_age + group_id, design.df = alldata_design)

da2_results %>%
  arrange(SpatialFDR) %>%
  head() 
##        logFC   logCPM        F       PValue       FDR Nhood SpatialFDR
##  1  6.326216 8.557685 18.46120 0.0001138818 0.2029138  1096  0.2205914
##  2  5.275664 8.200486 17.83941 0.0001440126 0.2029138  2423  0.2205914
##  3  4.837525 8.348456 15.57169 0.0002526991 0.2373687  2437  0.2629561
##  4  4.970073 8.332697 12.50975 0.0009028742 0.5749959   800  0.5956998
##  5 -5.872795 8.972963 12.32640 0.0010202198 0.5749959   819  0.5956998
##  6 -5.048635 9.391392 10.05721 0.0026315963 0.7318605  1038  0.7167789

p <- ggplot(da2_results, aes(PValue)) + geom_histogram(bins=50)
ggsave(file="mod/PValue_histo.DA2.pdf", plot=p, height=7, width=7)



da3_results <- testNhoods(alldata_milo, design = ~ group_id, design.df = alldata_design)

da3_results %>%
  arrange(SpatialFDR) %>%
  head() 
##       logFC   logCPM        F       PValue        FDR Nhood SpatialFDR
##  1 5.558100 8.348490 23.54180 3.026680e-05 0.08529183  2437 0.09822251
##  2 5.465724 8.626231 21.06748 6.466529e-05 0.09111340    25 0.10214123
##  3 6.822813 8.557656 22.00849 1.100066e-04 0.10333286  1096 0.11338843
##  4 6.152205 8.200521 20.61052 1.593314e-04 0.11224894  2423 0.12390216
##  5 4.452910 8.999469 15.29934 4.466782e-04 0.20978985  1796 0.20749301
##  6 4.439970 8.973483 15.48851 4.175209e-04 0.20978985  1968 0.20749301

p <- ggplot(da3_results, aes(PValue)) + geom_histogram(bins=50)
ggsave(file="mod/PValue_histo.DA3.pdf", plot=p, height=7, width=7)





#Volcano
p <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
ggsave(file="mod/DA_Volcano.pdf", plot=p, height=7, width=7)


p <- ggplot(da2_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
ggsave(file="mod/DA2_Volcano.pdf", plot=p, height=7, width=7)

p <- ggplot(da3_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
ggsave(file="mod/DA3_Volcano.pdf", plot=p, height=7, width=7)



alldata_milo <- buildNhoodGraph(alldata_milo)
## Plot single-cell UMAP
UMAP_pl <- plotReducedDim(alldata_milo, dimred = "UMAP", colour_by="group_id", text_by = "predicted.celltype.l2", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(alldata_milo, da_results, layout="UMAP",alpha=0.1)
nh2_graph_pl <- plotNhoodGraphDA(alldata_milo, da2_results, layout="UMAP",alpha=0.1)
nh3_graph_pl <- plotNhoodGraphDA(alldata_milo, da3_results, layout="UMAP",alpha=0.1)
##  Scale for fill is already present.
#@  Adding another scale for fill, which will replace the existing scale.

p <- UMAP_pl + nh_graph_pl
ggsave(file="mod/UMAP_DA_alpha0.1.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh2_graph_pl
ggsave(file="mod/UMAP_DA2_alpha0.1.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh3_graph_pl
ggsave(file="mod/UMAP_DA3_alpha0.1.png", plot=p, width=10, height=5)


# FDR 1
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, da_results, layout="UMAP",alpha=1)
nh2_graph_pl1 <- plotNhoodGraphDA(alldata_milo, da2_results, layout="UMAP",alpha=1)
nh3_graph_pl1 <- plotNhoodGraphDA(alldata_milo, da3_results, layout="UMAP",alpha=1)

p <- UMAP_pl + nh_graph_pl1
ggsave(file="mod/UMAP_DA_alpha1.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh2_graph_pl1
ggsave(file="mod/UMAP_DA2_alpha1.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh3_graph_pl1
ggsave(file="mod/UMAP_DA3_alpha1.png", plot=p, width=10, height=5)
ggsave(file="mod/UMAP_DA3_alpha1.pdf", plot=p, width=10, height=5)


# FDR 0.2
nh_graph_pl02 <- plotNhoodGraphDA(alldata_milo, da_results, layout="UMAP",alpha=0.2)
nh2_graph_pl02 <- plotNhoodGraphDA(alldata_milo, da2_results, layout="UMAP",alpha=0.2)
nh3_graph_pl02 <- plotNhoodGraphDA(alldata_milo, da3_results, layout="UMAP",alpha=0.2)

p <- UMAP_pl + nh_graph_pl02
ggsave(file="mod/UMAP_DA_alpha0.2.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh2_graph_pl02
ggsave(file="mod/UMAP_DA2_alpha0.2.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh3_graph_pl02
ggsave(file="mod/UMAP_DA3_alpha0.2.png", plot=p, width=10, height=5)

# FDR 0.4
nh_graph_pl04 <- plotNhoodGraphDA(alldata_milo, da_results, layout="UMAP",alpha=0.4)
nh2_graph_pl04 <- plotNhoodGraphDA(alldata_milo, da2_results, layout="UMAP",alpha=0.4)
nh3_graph_pl04 <- plotNhoodGraphDA(alldata_milo, da3_results, layout="UMAP",alpha=0.4)

p <- UMAP_pl + nh_graph_pl04
ggsave(file="mod/UMAP_DA_alpha0.4.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh2_graph_pl04
ggsave(file="mod/UMAP_DA2_alpha0.4.png", plot=p, width=10, height=5)

p <- UMAP_pl + nh3_graph_pl04
ggsave(file="mod/UMAP_DA3_alpha0.4.png", plot=p, width=10, height=5)


### annotate
da_results <- annotateNhoods(alldata_milo, da_results, coldata_col = "predicted.celltype.l2")
da2_results <- annotateNhoods(alldata_milo, da2_results, coldata_col = "predicted.celltype.l2")
da3_results <- annotateNhoods(alldata_milo, da3_results, coldata_col = "predicted.celltype.l2")

da_results <- annotateNhoods(alldata_milo, da_results, coldata_col = "SCT_snn_res.0.4")
da2_results <- annotateNhoods(alldata_milo, da2_results, coldata_col = "SCT_snn_res.0.4")
da3_results <- annotateNhoods(alldata_milo, da3_results, coldata_col = "SCT_snn_res.0.4")


da3_results %>%
  arrange(SpatialFDR) %>%
  head() 
##       logFC   logCPM        F       PValue        FDR Nhood SpatialFDR
##  1 5.558100 8.348490 23.54180 3.026680e-05 0.08529183  2437 0.09822251
##  2 5.465724 8.626231 21.06748 6.466529e-05 0.09111340    25 0.10214123
##  3 6.822813 8.557656 22.00849 1.100066e-04 0.10333286  1096 0.11338843
##  4 6.152205 8.200521 20.61052 1.593314e-04 0.11224894  2423 0.12390216
##  5 4.452910 8.999469 15.29934 4.466782e-04 0.20978985  1796 0.20749301
##  6 4.439970 8.973483 15.48851 4.175209e-04 0.20978985  1968 0.20749301
##    predicted.celltype.l2 predicted.celltype.l2_fraction SCT_snn_res.0.4
##  1               CD8 TEM                       1.000000              17
##  2               CD8 TEM                       1.000000              17
##  3             CD14 Mono                       1.000000               1
##  4               CD8 TEM                       1.000000              17
##  5             CD14 Mono                       1.000000               1
##  6             CD14 Mono                       0.984127               1
##    SCT_snn_res.0.4_fraction
##  1                0.9701493
##  2                0.9784946
##  3                1.0000000
##  4                1.0000000
##  5                1.0000000
##  6                1.0000000



p <- ggplot(da3_results, aes(predicted.celltype.l2_fraction)) + geom_histogram(bins=50) 
ggsave(file="mod/DA3_celltype_fraction.pdf", plot=p, width=7, height=7)

p <- ggplot(da3_results, aes(SCT_snn_res.0.4_fraction)) + geom_histogram(bins=50) 
ggsave(file="mod/DA3_SCT_snn_res.0.4_fraction.pdf", plot=p, width=7, height=7)

da3_results$predicted.celltype.l2 <- ifelse(da3_results$predicted.celltype.l2_fraction < 0.8, "Mixed", da3_results$predicted.celltype.l2)

p <- plotDAbeeswarm(da3_results, group.by = "predicted.celltype.l2", alpha = 1)
ggsave(file="mod/DA3_lFC_celltype.Mix0.8.pdf", plot=p, width=7, height=7)
##  Converting group.by to factor...
##  Scale for colour is already present.
##  Adding another scale for colour, which will replace the existing scale.
##  Warning message:
##  Removed 2899 rows containing missing values.

p <- plotDAbeeswarm(da3_results, group.by = "predicted.celltype.l2", alpha = 0.1)
ggsave(file="mod/DA3_lFC_celltypeMix0.8_alpha0.1.pdf", plot=p, width=7, height=7)





### Finding markers of DA populations
## Automatic grouping of neighbourhoods
## Run buildNhoodGraph to store nhood adjacency matrix
alldata_milo <- buildNhoodGraph(alldata_milo)

## Find groups
da3_results <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 2)
##  Found 1 DA neighbourhoods at FDR 10%
##  nhoodAdjacency found - using for nhood grouping
head(da3_results)
##         logFC   logCPM          F    PValue       FDR Nhood SpatialFDR
##  1  0.1585037 8.241322 0.02692985 0.8706622 0.9999983     1  0.9999983
##  2  0.5801598 8.660933 0.32826609 0.5706704 0.9999983     2  0.9868911
##  3 -0.7416344 8.704363 0.51328071 0.4788979 0.9999983     3  0.9868911
##  4  0.1886788 8.634268 0.04225085 0.8384194 0.9999983     4  0.9999983
##  5 -1.6766969 8.645994 1.78444228 0.1909909 0.9880502     5  0.9476888
##  6 -0.4287207 9.118397 0.15895046 0.6927621 0.9999983     6  0.9969338
##    predicted.celltype.l2 predicted.celltype.l2_fraction SCT_snn_res.0.4
##  1               CD4 TCM                      0.9250000               5
##  2                 Mixed                      0.6986301               0
##  3                 Mixed                      0.7812500               7
##  4                 Mixed                      0.4918033               6
##  5               CD8 TEM                      0.8135593               0
##  6                 Mixed                      0.7676768               2
##    SCT_snn_res.0.4_fraction NhoodGroup
##  1                 0.775000          1
##  2                 1.000000          2
##  3                 0.781250          3
##  4                 0.852459          1
##  5                 1.000000          2
##  6                 1.000000          4

p <- plotNhoodGroups(alldata_milo, da3_results, layout="UMAP") 
ggsave(file="mod/DA3_NhoodGroups.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(da3_results, group.by = "NhoodGroup", alpha = 1)
ggsave(file="mod/DA3_lFC_NhoodGroups.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 0.5), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=0.5")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta0.5.pdf", plot=p, width=7, height=7)

set.seed(42)
p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 1), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=1")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta1.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 1.5), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=1.5")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta1.5.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 2), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=2")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta2.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 3), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=3")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta3.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 4), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=4")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta4.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 5), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=5")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta5.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 3, overlap=1), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=3 & overlap=1")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta3.overlap1.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 3, overlap=3), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=3 & overlap=3")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta3.overlap3.pdf", plot=p, width=7, height=7)

p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 3, overlap=5), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=3 & overlap=5")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta3.overlap5.pdf", plot=p, width=7, height=7)


# LFC=100
set.seed(42)
p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=1), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=100 & overlap=1")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta100.overlap1.pdf", plot=p, width=7, height=7)
set.seed(42)
p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=5), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=100 & overlap=5")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta100.overlap5.pdf", plot=p, width=7, height=7)
set.seed(42)
p <- plotDAbeeswarm(groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=6), group.by = "NhoodGroup", alpha=1) + ggtitle("max LFC delta=100 & overlap=6")
ggsave(file="mod/DA3_lFC_NhoodGroups.delta100.overlap6.pdf", plot=p, width=7, height=7)



### Let’s settle for overlap=1 and max.lfc.delta=3
set.seed(42)
da3_results <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 3, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta3_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta3_overlap1.png", plot=p, width=7, height=7)


set.seed(42)
da3_results_delta2 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 2, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results_delta2, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta2_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta2_overlap1.png", plot=p, width=7, height=7)

set.seed(42)
da3_results_delta1 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 1, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results_delta1, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta1_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta1_overlap1.png", plot=p, width=7, height=7)

set.seed(42)
da3_results_delta1.5 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 1.5, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results_delta1.5, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta1.5_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta1.5_overlap1.png", plot=p, width=7, height=7)

set.seed(42)
da3_results_delta0.5 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 0.5, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results_delta0.5, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta0.5_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta0.5_overlap1.png", plot=p, width=7, height=7)

set.seed(42)
da3_results_delta5 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 5)
p <- plotNhoodGroups(alldata_milo, da3_results_delta5, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta5_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta5_overlap1.png", plot=p, width=7, height=7)


set.seed(42)
da3_results_overlap1 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=1)
p <- plotNhoodGroups(alldata_milo, da3_results_overlap1, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap1.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap1.png", plot=p, width=7, height=7)

set.seed(42)
da3_results_overlap5 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=5)
p <- plotNhoodGroups(alldata_milo, da3_results_overlap5, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap5.pdf", plot=p, width=6.5, height=5)
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap5.png", plot=p, width=6.5, height=5)

set.seed(42)
da3_results_overlap6 <- groupNhoods(alldata_milo, da3_results, max.lfc.delta = 100, overlap=6)
p <- plotNhoodGroups(alldata_milo, da3_results_overlap6, layout="UMAP")
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap6.pdf", plot=p, width=7, height=7)
ggsave(file="mod/DA3_NhoodGroups_delta100_overlap6.png", plot=p, width=7, height=7)








### GO TO milo_mod_overlap5.R

######################################
