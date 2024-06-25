seurat.combined.miQCed.Scdsed.filt <- readRDS("/mnt/e/scRNA_SLE/seurat_qc/seurat.combined.miQCed.Scdsed.filt")

library(harmony)
library(Seurat)


### SplitObject
seurat.list <- SplitObject(object = seurat.combined.miQCed.Scdsed.filt, split.by = "Sample")

### SCTransform normalization separately
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- SCTransform(seurat.list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = TRUE)
}

### Integration Features without V/D/J genes & no_named_gene       # 予定のn_featuresより1000多くしている(V/D/J遺伝子などを除去するため)
for (i in 1:length(seurat.list)) {
    DefaultAssay(seurat.list[[i]]) <- "SCT"
}
seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)

library(tidyverse)
seurat.features_tibble <- seurat.features %>% as_tibble
out_gene <- c("^TRAV","^TRAJ","^TRBV","^TRDV","^TRBJ","^IGH","^IGL","^IGK")                                              # B cell immunoglobulin data: IGH, IGK, IGL
seurat.features_tibble %>% filter(nchar(value) == 10) %>% filter(substr(.$value,9,9) == ".") %>% .$value -> noname_gene   # no_named_genes 
seurat.features_tibble %>% filter(!str_detect(value, paste(out_gene, collapse = "|"))) %>% filter(!str_detect(value, paste(noname_gene, collapse = "|"))) %>% .$value -> gene_filtered

seurat.features <- gene_filtered[1:2000]

write(seurat.features, file = "2000_features_fromSCT.txt") 

saveRDS(seurat.list, "seurat.list_SCTed")
# seurat.list <- readRDS("seurat.list_SCTed")



### Azimuth annotation per Sample
library(Azimuth)

# Load the reference
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "/mnt/e/scRNA_SLE/azimuth_references/human_pbmc/v1.0.0")

# Preprocess with SCTransform
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- SCTransform(
        object = seurat.list[[i]],
        assay = "RNA",
        new.assay.name = "refAssay",
        residual.features = rownames(x = reference$map),
        reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
        method = 'glmGamPoi',
        ncells = 2000,
        n_genes = 2000,
        do.correct.umi = FALSE,
        do.scale = FALSE,
        do.center = TRUE
        )    
}

# Find anchors between query and reference
anchors <- list()
for (i in 1:length(seurat.list)) {
    anchors[[i]] <- FindTransferAnchors(
        reference = reference$map,
        query = seurat.list[[i]],
        k.filter = NA,
        reference.neighbors = "refdr.annoy.neighbors",
        reference.assay = "refAssay",
        query.assay = "refAssay",
        reference.reduction = "refDR",
        normalization.method = "SCT",
        features = intersect(rownames(x = reference$map), VariableFeatures(object = seurat.list[[i]])),
        dims = 1:50,
        n.trees = 20,
        mapping.score.k = 100
        )
}

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = c("celltype.l2", "celltype.l1", "celltype.l3"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("celltype.l2", "celltype.l1", "celltype.l3")
if (TRUE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}

for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- TransferData(
        reference = reference$map,
        query = seurat.list[[i]],
        dims = 1:50,
        anchorset = anchors[[i]],
        refdata = refdata,
        n.trees = 20,
        store.weights = TRUE
        )
}

# Calculate the embeddings of the query data on the reference SPCA
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- IntegrateEmbeddings(
        anchorset = anchors[[i]],
        reference = reference$map,
        query = seurat.list[[i]],
        reductions = "pcaproject",
        reuse.weights.matrix = TRUE
        )
}

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
for (i in 1:length(seurat.list)) {
    seurat.list[[i]][["query_ref.nn"]] <- FindNeighbors(
        object = Embeddings(reference$map[["refDR"]]),
        query = Embeddings(seurat.list[[i]][["integrated_dr"]]),
        return.neighbor = TRUE,
        l2.norm = TRUE
        )
}

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.

#' Transform an NN index
#'
#' @param object Seurat object
#' @param meta.data Metadata
#' @param neighbor.slot Name of Neighbor slot
#' @param key Column of metadata to use
#'
#' @return \code{object} with transfomed neighbor.slot
#'
#' @importFrom SeuratObject Indices
#'
#' @keywords internal
#'
NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- SeuratObject::Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- NNTransform(
        object = seurat.list[[i]],
        meta.data = reference$map[[]]
        )
}

# Project the query to the reference UMAP.
for (i in 1:length(seurat.list)) {
    seurat.list[[i]][["proj.umap"]] <- RunUMAP(
        object = seurat.list[[i]][["query_ref.nn"]],
        reduction.model = reference$map[["refUMAP"]],
        reduction.key = 'UMAP_'
        )
}

# Calculate mapping score and add to metadata
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- AddMetaData(
        object = seurat.list[[i]],
        metadata = MappingScore(anchors = anchors[[i]]),
        col.name = "mapping.score"
        )
}

# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
celltypes <- c("celltype.l2", "celltype.l1", "celltype.l3")

dir.create("Azimuth_Plot", recursive = T)
dir.create("Azimuth_Plot/Prediction_Score", recursive = T)
dir.create("Azimuth_Plot/Mapping_Score", recursive = T)
dir.create("Azimuth_Plot/DimPlot", recursive = T)
dir.create("Azimuth_Plot/RNA", recursive = T)


for (n in 1:length(celltypes)) {
    id <- celltypes[[n]]
    predicted.id <- paste0("predicted.", id)

    # DimPlot of the reference
    png(paste0("Azimuth_Plot/DimPlot/reference.",id,".png"), width = 800, height = 500)
    p <- DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE, repel=TRUE)
    p
    dev.off()

    for (i in 1:length(seurat.list)) {
        sampleid <- names(seurat.list[i])

        # DimPlot of the query, colored by predicted cell type
        png(paste0("Azimuth_Plot/DimPlot/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
        p <- DimPlot(object = seurat.list[[i]], reduction = "proj.umap", group.by = predicted.id, label = TRUE, repel=TRUE)
        plot(p)
        dev.off()
        
        # Plot the score for the predicted cell type of the query
        png(paste0("Azimuth_Plot/Prediction_Score/",sampleid,".",predicted.id,".png"), width = 1000, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = paste0(predicted.id, ".score"), reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()
        plot(p)
        dev.off()
        
        # Plot the mapping score
        png(paste0("Azimuth_Plot/Mapping_Score/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = "mapping.score", reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = "mapping.score", group.by = predicted.id) + NoLegend()
        plot(p)
        dev.off()
        
        ## Plot the prediction score for the class CD16 Mono
        #dir.create(paste0("Azimuth_Plot/Prediction_Score/","CD16_Mono"), recursive = T)
        #png(paste0("Azimuth_Plot/Prediction_Score/","CD16_Mono","/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
        #FeaturePlot(object = seurat.list[[i]], features = "CD16 Mono", reduction = "proj.umap")
        #VlnPlot(object = seurat.list[[i]], features = "CD16 Mono", group.by = predicted.id) + NoLegend()
        #dev.off()
        
        # Plot an RNA feature
        dir.create(paste0("Azimuth_Plot/RNA/","GNLY"), recursive = T)
        png(paste0("Azimuth_Plot/RNA/","GNLY","/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = "GNLY", reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = "GNLY", group.by = predicted.id) + NoLegend()
        dev.off()
        
        # Plot an imputed protein feature
        if (TRUE) {
            dir.create(paste0("Azimuth_Plot/Imputed_Protein/","CD3-1"), recursive = T)
            png(paste0("Azimuth_Plot/Imputed_Protein/","CD3-1","/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
            p <- FeaturePlot(object = seurat.list[[i]], features = "CD3-1", reduction = "proj.umap")
            p <- p + VlnPlot(object = seurat.list[[i]], features = "CD3-1", group.by = predicted.id) + NoLegend()
            plot(p)
            dev.off()
            }
        }
}



### merge
seurat_merged <- merge(x = seurat.list[[1]], y = seurat.list[2:length(seurat.list)], project = "SLE", merge.data = TRUE)
DefaultAssay(seurat_merged) <- "SCT"


### Manually set variable features of merged Seurat object        highly variable gene (HVG)
VariableFeatures(seurat_merged) <- seurat.features

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_merged), 10)

vdj_gene <- seurat.features_tibble %>% filter(str_detect(value, paste(out_gene, collapse = "|"))) %>% .$value
repel_gene = c(top10, vdj_gene, noname_gene)

# plot variable features with and without labels
for (i in 1:length(seurat.list)) {
    plot1 <- VariableFeaturePlot(seurat.list[[i]], assay = "SCT", selection.method="sct")
    plot_final <- paste0("plot_", i)
    assign(plot_final, LabelPoints(plot = plot1, points = repel_gene, repel = TRUE, xnudge = 0, ynudge = 0))
}
png("SCT.VariableFeatures.png", width=2000, height=2000)
plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 + plot_9 + plot_10 + plot_11 + plot_12 + plot_13 + plot_14 + plot_15
dev.off()



### PCA
seurat_merged <- RunPCA(object = seurat_merged, assay = "SCT", npcs = 50)


png("SCT.VizDimLoadings.png", width = 800, height = 500)
VizDimLoadings(seurat_merged, dims = 1:2, reduction = "pca")
dev.off()

png("SCT.pca.Sample.png", width = 800, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="Sample")
dev.off()

png("SCT.pca.Institute.png", width = 800, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="Institute")
dev.off()

png("SCT.pca.eHHV-6B.png", width = 800, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="SampleGroup")
dev.off()

png("SCT.pca.png", width = 800, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="Institute", split.by="SampleGroup")
dev.off()

png("SCT.heatmap_dim_1_15.png", width = 800, height = 1500)
DimHeatmap(seurat_merged, dims = 1:15, cells = 5000, balanced = TRUE)
dev.off()

### SCTでは各遺伝子の分散が等しいという仮定が上手くいかないのでJackStrawは動作しないようになっている。　https://github.com/satijalab/seurat/issues/1871

### ElbowPlot
png("SCT.ElbowPlot.png", width=800, height=500)
ElbowPlot(seurat_merged, ndims = 50, reduction = "pca")
dev.off()


### Harmony         ElbowPlotをみてとりあえずPC50まで全部にした
seurat_merged <- RunHarmony(object = seurat_merged,
                                    assay.use = "SCT",
                                    reduction = "pca",
                                    dims.use = 1:50,
                                    group.by.vars = c("Institute", "Sample"),
                                    plot_convergence = TRUE)
##  Harmony converged after 8 iterations
##  Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity


### UMAP with harmony       
seurat_merged <- RunUMAP(object = seurat_merged, assay = "SCT", reduction = "harmony", dims = 1:30)
seurat_merged <- FindNeighbors(object = seurat_merged, assay = "SCT", reduction = "harmony", dims = 1:30)
seurat_merged <- FindClusters(object = seurat_merged, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.2))


### Change default assay to "RNA"; normalize then generate FeaturePlots and perform differential expression analysis    https://github.com/satijalab/seurat/issues/1836
DefaultAssay(seurat_merged) <-"RNA"
seurat_merged <- NormalizeData(seurat_merged, verbose = TRUE) 

saveRDS(seurat_merged,"Harmony_round1.rds")

seurat_merged
##  An object of class Seurat
##  76712 features across 66915 samples within 7 assays
##  Active assay: RNA (27204 features, 0 variable features)
##   6 other assays present: SCT, refAssay, prediction.score.celltype.l2, prediction.score.celltype.l1, prediction.score.celltype.l3, impADT
##   3 dimensional reductions calculated: pca, harmony, umap

seurat_merged[["refAssay"]] <- NULL
seurat_merged[["prediction.score.celltype.l2"]] <- NULL
seurat_merged[["prediction.score.celltype.l1"]] <- NULL
seurat_merged[["prediction.score.celltype.l3"]] <- NULL
seurat_merged[["impADT"]] <- NULL

seurat_merged
##  An object of class Seurat
##  49185 features across 66915 samples within 2 assays
##  Active assay: RNA (27204 features, 0 variable features)
##   1 other assay present: SCT
##   3 dimensional reductions calculated: pca, harmony, umap

library(SeuratDisk)
SaveH5Seurat(seurat_merged, filename = "Harmony_round1.h5Seurat")

Convert("Harmony_round1.h5Seurat", dest = "h5ad")

seurat_merged[["SCT"]] <- NULL

seurat_merged
##  An object of class Seurat
##  27204 features across 66915 samples within 1 assay
##  Active assay: RNA (27204 features, 0 variable features)
##   1 dimensional reduction calculated: umap

saveRDS(seurat_merged,"Harmony_RNA_round1.rds")

library(SeuratDisk)
SaveH5Seurat(seurat_merged, filename = "Harmony_RNA_round1.h5Seurat")

Convert("Harmony_RNA_round1.h5Seurat", dest = "h5ad")



#-------------------------------#
#             graph             #
#-------------------------------#            
p_list <- list()
for(meta in c("Sample","SampleGroup","Institute")){
    seurat_merged <- SetIdent(seurat_merged, value = meta)
    p <- DimPlot(seurat_merged, reduction = "umap", group.by = meta) +
            theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
    ggsave(file = paste0("UMAP_", meta,".png"), plot = p, width =12, height = 10, dpi=220)
    
    p_list[[meta]] <- p

    }

for(meta in c("SCT_snn_res.0.2","SCT_snn_res.0.4","predicted.celltype.l1","predicted.celltype.l2","predicted.celltype.l3")){
    seurat_merged <- SetIdent(seurat_merged, value = meta)
    p <- DimPlot(seurat_merged, reduction = "umap", group.by = meta) + NoLegend() +
             theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
    p <- LabelClusters(plot = p, id = meta, size=4.5)
    ggsave(file = paste0("UMAP_", meta,".png"), plot = p, width =10, height = 10, dpi=220)
    
    p_list[[meta]] <- p

    }


p_seurat_merged <- (p_list[["predicted.celltype.l1"]] + p_list[["predicted.celltype.l2"]]) /
           (p_list[["SCT_snn_res.0.2"]] + p_list[["SCT_snn_res.0.4"]])

ggsave(file = "UMAP_clustering_seurat_merged.png", plot = p_seurat_merged, width =15, height = 12.5, dpi=220)


# FeaturePlot
fp.plt <- FeaturePlot(seurat_merged, reduction = "umap", c("PPBP", "GP9"))
ggsave(file = "FeaturePlot_plt.png", plot = fp.plt , width = 10, height = 4.5, dpi = 220)

fp.Innate <- FeaturePlot(seurat_merged, reduction = "umap", c("CD14", "FCGR3A","CD68", "LYZ"))
ggsave(file = "FeaturePlot_Innate.png", plot = fp.Innate , width = 10, height = 9, dpi = 220)

fp.DC <- FeaturePlot(seurat_merged, reduction = "umap", c("CD1C","CLEC10A","LILRA4","IL3RA"))
ggsave(file = "FeaturePlot_DC.png", plot = fp.DC , width = 10, height = 9, dpi = 220)

fp.T_NK <- FeaturePlot(seurat_merged, reduction = "umap", c("CD3E", "CD8A", "CD4", "KLRF1"))
ggsave(file = "FeaturePlot_T_NK.png", plot = fp.T_NK , width = 10, height = 9, dpi = 220)

fp.T2 <- FeaturePlot(seurat_merged, reduction = "umap", c("FOXP3","SLC4A10","TRAV1-2", "MKI67"))
ggsave(file = "FeaturePlot_T_minor.png", plot = fp.T2 , width = 10, height = 9, dpi = 220)

fp.B <- FeaturePlot(seurat_merged, reduction = "umap", c("CD79A", "MS4A1", "MZB1", "S100A9"))           # S100A9: Doublets (Plasmablast & Monocytes) 
ggsave(file = "FeaturePlot_B.png", plot = fp.B , width = 10, height = 9, dpi = 220)



p <- DimPlot(seurat_merged, reduction="umap", split.by="Sample", group.by="predicted.celltype.l2", ncol=4) + NoLegend() +
             theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
p <- LabelClusters(plot = p, id = "predicted.celltype.l2", size=4.5)
ggsave(file="UMAP_perSample.png", plot=p, width=30, height=30, dpi=220, limitsize=FALSE)


### CD14 vs CD16
p <- FeaturePlot(seurat_merged, reduction = "umap", c("CD14", "FCGR3A"), blend=TRUE)
ggsave(file = "FeaturePlot_mono.png", plot = p , width = 30, height = 9, dpi = 220)


# TNFRSF4 (aka CD134 or OX40, the receptor for HHV-6B) (paper: Latent human herpesvirus 6 is reaactivated in CAR T cells)
seurat_merged <- SetIdent(seurat_merged, value = "predicted.celltype.l2")
p <- FeaturePlot(seurat_merged, reduction = "umap", features = "TNFRSF4", split.by="SampleGroup", label=TRUE, repel = TRUE, keep.scale="feature") &
        theme(legend.position = c(0.95, 0.1))
ggsave(file = "FeaturePlot_TNFRSF4.png", plot = p , width = 20, height = 9, dpi = 220)



### antiviral ISG score     from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3274382/  
# review https://www.annualreviews.org/doi/10.1146/annurev-virology-092818-015756?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed
# review https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6334645/
antiviral_ISGs <- list(c(
    'ISG15',
    'ISG20',
    'ADARB1',
    'APOBEC3A',
    'APOBEC3G',
    'APOBEC3H',
    'APOBEC3C',
    'CD74',
    'DDIT4',
    'DDX60',
    'GBP1',
    'GBP2',
    'GBP5',
    'GBP4',
    'IFI44L',
    'IFI6',
    'IFIT3',
    'IFIT2',
    'IFIT1',
    'IFITM3',
    'IFITM1',
    'IFITM2',
    'LY6E',
    'IFI27',
    'HERC5',
    'LGALS3BP',
    'OAS1',
    'OAS3',
    'OASL',
    'RIN2',
    'SAMHD1',
    'SERINC5',
    'RSAD2'
    ))
seurat_merged <- AddModuleScore(
    object = seurat_merged,
    features = antiviral_ISGs,
    ctrl = 100,
    name = 'antiviral_ISGs'
    )

library(RColorBrewer)
seurat_merged <- SetIdent(seurat_merged, value = "predicted.celltype.l2")
p <- FeaturePlot(seurat_merged, reduction = "umap", features = "antiviral_ISGs1", split.by="SampleGroup", label=TRUE, repel = TRUE, keep.scale="feature") &
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
        theme(legend.position = c(0.95, 0.1))
ggsave(file = "FeaturePlot_antiviralISGs.png", plot = p , width = 20, height = 9, dpi = 220)

### type-1 ISG signature    https://www.nature.com/articles/s41467-022-35209-1
ISGs <- list(c(
    'ISG15',
    'IFI6',
    'IFI44L',
    'IFI44',
    'RSAD2',
    'CXCL10',
    'IFIT2',
    'IFIT3',
    'IFIT1',
    'IFITM3',
    'OAS1',
    'OAS3',
    'OASL',
    'EPST11',
    'RNASE1',
    'RNASE2',
    'IFI27',
    'XAF1',
    'IGALS3BP',
    'SIGLEC1',
    'USP18',
    'APOBEC3A',
    'APOBEC3B',
    'MX1'
    ))
seurat_merged <- AddModuleScore(
    object = seurat_merged,
    features = ISGs,
    ctrl = 100,
    name = 'type-1_ISGs'
    )

library(RColorBrewer)
seurat_merged <- SetIdent(seurat_merged, value = "predicted.celltype.l2")
p <- FeaturePlot(seurat_merged, reduction = "umap", features = "type.1_ISGs1", split.by="SampleGroup", label=TRUE, repel = TRUE, keep.scale="feature") &
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
        theme(legend.position = c(0.95, 0.1))
ggsave(file = "FeaturePlot_type-1ISGs.png", plot = p , width = 20, height = 9, dpi = 220)



saveRDS(seurat_merged,"Harmony_round1.ISGs.rds")
