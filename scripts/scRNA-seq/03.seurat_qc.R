library(data.table)
library(Seurat)
library(miQC)
library(SeuratWrappers)
library(flexmix)
library(SingleCellExperiment)
library(Matrix)
library(stringr)


sce01 <- readRDS("../scds/count_id01/sce")
sce02 <- readRDS("../scds/count_id02/sce")
sce03 <- readRDS("../scds/count_id03/sce")
sce04 <- readRDS("../scds/count_id04/sce")
sce05 <- readRDS("../scds/count_id05/sce")
sce06 <- readRDS("../scds/count_id06/sce")
sce07 <- readRDS("../scds/count_id07/sce")
sce08 <- readRDS("../scds/count_id08/sce")
sce09 <- readRDS("../scds/count_id09/sce")
sce10 <- readRDS("../scds/count_id10/sce")
sce11 <- readRDS("../scds/count_id11/sce")
sce12 <- readRDS("../scds/count_id12/sce")
sce13 <- readRDS("../scds/count_id13/sce")
sce14 <- readRDS("../scds/count_id14/sce")
sce15 <- readRDS("../scds/count_id15/sce")



### Create SeuratObject
seurat01 <- CreateSeuratObject(counts=counts(sce01))
seurat02 <- CreateSeuratObject(counts=counts(sce02))
seurat03 <- CreateSeuratObject(counts=counts(sce03))
seurat04 <- CreateSeuratObject(counts=counts(sce04))
seurat05 <- CreateSeuratObject(counts=counts(sce05))
seurat06 <- CreateSeuratObject(counts=counts(sce06))
seurat07 <- CreateSeuratObject(counts=counts(sce07))
seurat08 <- CreateSeuratObject(counts=counts(sce08))
seurat09 <- CreateSeuratObject(counts=counts(sce09))
seurat10 <- CreateSeuratObject(counts=counts(sce10))
seurat11 <- CreateSeuratObject(counts=counts(sce11))
seurat12 <- CreateSeuratObject(counts=counts(sce12))
seurat13 <- CreateSeuratObject(counts=counts(sce13))
seurat14 <- CreateSeuratObject(counts=counts(sce14))
seurat15 <- CreateSeuratObject(counts=counts(sce15))

### add metadata
seurat01$Sample <- sce01$Sample
seurat02$Sample <- sce02$Sample
seurat03$Sample <- sce03$Sample
seurat04$Sample <- sce04$Sample
seurat05$Sample <- sce05$Sample
seurat06$Sample <- sce06$Sample
seurat07$Sample <- sce07$Sample
seurat08$Sample <- sce08$Sample
seurat09$Sample <- sce09$Sample
seurat10$Sample <- sce10$Sample
seurat11$Sample <- sce11$Sample
seurat12$Sample <- sce12$Sample
seurat13$Sample <- sce13$Sample
seurat14$Sample <- sce14$Sample
seurat15$Sample <- sce15$Sample

seurat01$Barcode <- sce01$Barcode
seurat02$Barcode <- sce02$Barcode
seurat03$Barcode <- sce03$Barcode
seurat04$Barcode <- sce04$Barcode
seurat05$Barcode <- sce05$Barcode
seurat06$Barcode <- sce06$Barcode
seurat07$Barcode <- sce07$Barcode
seurat08$Barcode <- sce08$Barcode
seurat09$Barcode <- sce09$Barcode
seurat10$Barcode <- sce10$Barcode
seurat11$Barcode <- sce11$Barcode
seurat12$Barcode <- sce12$Barcode
seurat13$Barcode <- sce13$Barcode
seurat14$Barcode <- sce14$Barcode
seurat15$Barcode <- sce15$Barcode

seurat01$SampleGroup <- sce01$SampleGroup
seurat02$SampleGroup <- sce02$SampleGroup
seurat03$SampleGroup <- sce03$SampleGroup
seurat04$SampleGroup <- sce04$SampleGroup
seurat05$SampleGroup <- sce05$SampleGroup
seurat06$SampleGroup <- sce06$SampleGroup
seurat07$SampleGroup <- sce07$SampleGroup
seurat08$SampleGroup <- sce08$SampleGroup
seurat09$SampleGroup <- sce09$SampleGroup
seurat10$SampleGroup <- sce10$SampleGroup
seurat11$SampleGroup <- sce11$SampleGroup
seurat12$SampleGroup <- sce12$SampleGroup
seurat13$SampleGroup <- sce13$SampleGroup
seurat14$SampleGroup <- sce14$SampleGroup
seurat15$SampleGroup <- sce15$SampleGroup

seurat01$Institute <- sce01$Institute
seurat02$Institute <- sce02$Institute
seurat03$Institute <- sce03$Institute
seurat04$Institute <- sce04$Institute
seurat05$Institute <- sce05$Institute
seurat06$Institute <- sce06$Institute
seurat07$Institute <- sce07$Institute
seurat08$Institute <- sce08$Institute
seurat09$Institute <- sce09$Institute
seurat10$Institute <- sce10$Institute
seurat11$Institute <- sce11$Institute
seurat12$Institute <- sce12$Institute
seurat13$Institute <- sce13$Institute
seurat14$Institute <- sce14$Institute
seurat15$Institute <- sce15$Institute

seurat01$Scds_call <- sce01$hybrid_call
seurat02$Scds_call <- sce02$hybrid_call
seurat03$Scds_call <- sce03$hybrid_call
seurat04$Scds_call <- sce04$hybrid_call
seurat05$Scds_call <- sce05$hybrid_call
seurat06$Scds_call <- sce06$hybrid_call
seurat07$Scds_call <- sce07$hybrid_call
seurat08$Scds_call <- sce08$hybrid_call
seurat09$Scds_call <- sce09$hybrid_call
seurat10$Scds_call <- sce10$hybrid_call
seurat11$Scds_call <- sce11$hybrid_call
seurat12$Scds_call <- sce12$hybrid_call
seurat13$Scds_call <- sce13$hybrid_call
seurat14$Scds_call <- sce14$hybrid_call
seurat15$Scds_call <- sce15$hybrid_call




### merge (not integrate)
seurat.combined <- merge(seurat01, y = c(seurat02,seurat03,seurat04,seurat05,seurat06,seurat07,seurat08,seurat09,seurat10,seurat11,seurat12,seurat13,seurat14,seurat15,seurat16), add.cell.ids = c(sce01$Sample[1], sce02$Sample[1], sce03$Sample[1], sce04$Sample[1], sce05$Sample[1], sce06$Sample[1], sce07$Sample[1], sce08$Sample[1], sce09$Sample[1], sce10$Sample[1], sce11$Sample[1], sce12$Sample[1], sce13$Sample[1], sce14$Sample[1], sce15$Sample[1]), project = "n15")

### QC
# MT
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^MT-", col.name = "percent.mt")
# RIBO
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^RP[SL]", col.name = "percent.ribo")
# Hb
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^HB[^(P)]", col.name = "percent.hb")
# Plat
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "PECAM1|PF4", col.name = "percent.plat")


png("QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0.1)
dev.off()

png("QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0.1)
dev.off()

png("QC_metrics.violin.Institute.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = TRUE)
dev.off()

png("QC_metrics.violin.Institute.Doublet.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = TRUE)
dev.off()

png("QC_metrics.violin.Institute.Doublet.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = FALSE)
dev.off()

png("QC_metrics.violin.Institute.miQC.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="miQC.keep", log = TRUE)
dev.off()

png("QC_metrics.violin.Institute.miQC.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="miQC.keep", log = FALSE)
dev.off()

png("QC_metrics.scatter.nCount_RNA.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

png("QC_metrics.scatter.nFeature_RNA.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

png("QC_metrics.scatter.percent.mt.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

png("QC_metrics.scatter.percent.ribo.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

png("QC_metrics.scatter.percent.hb.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

png("QC_metrics.scatter.percent.plat.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5
dev.off()

# filter out low quality cells using miQC
seurat.combined <- RunMiQC(seurat.combined, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, model.slot = "flexmix_model")
#   Warning: Adding a command log without an assay associated with it

pdf("PlotMiQC.pdf")
p <- PlotMiQC(seurat.combined, percent.mt="percent.mt", nFeature_RNA="nFeature_RNA", model.slot="flexmix_model", color.by="miQC.probability")
p
dev.off()

seurat.combined.miQCed <- subset(seurat.combined, miQC.keep == "keep")

png("miQCed.QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0.1)
dev.off()

png("miQCed.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0.1)
dev.off()

png("miQCed.QC_metrics.violin.Institute.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = TRUE)
dev.off()

png("miQCed.QC_metrics.violin.Institute.Doublet.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = TRUE)
dev.off()

png("miQCed.QC_metrics.violin.Institute.Doublet.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = FALSE)
dev.off()



seurat.combined.miQCed.Scdsed <- subset(seurat.combined.miQCed, Scds_call == "FALSE")

png("miQCed.Scdsed.QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0.1)
dev.off()

png("miQCed.Scdsed.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0.1)
dev.off()

png("miQCed.Scdsed.QC_metrics.violin.Institute.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = TRUE)
dev.off()

png("miQCed.Scdsed.QC_metrics.violin.Institute.Doublet.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = TRUE)
dev.off()

png("miQCed.Scdsed.QC_metrics.violin.Institute.Doublet.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = FALSE)
dev.off()


selected_gene <- rownames(seurat.combined.miQCed.Scdsed)[Matrix::rowSums(seurat.combined.miQCed.Scdsed) > 3]

seurat.combined.miQCed.Scdsed.filt <- subset(seurat.combined.miQCed.Scdsed, features=selected_gene, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 60000 & percent.ribo > 5 & percent.hb < 10)

png("miQCed.Scdsed.filt.QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0.1)
dev.off()

png("miQCed.Scdsed.filt.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0.1)
dev.off()

png("miQCed.Scdsed.filt.QC_metrics.violin.Institute.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = TRUE)
dev.off()

png("miQCed.Scdsed.filt.QC_metrics.violin.Institute.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined.miQCed.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = FALSE)
dev.off()


saveRDS(seurat.combined, file="seurat.combined")
saveRDS(seurat.combined.miQCed.Scdsed.filt, file="seurat.combined.miQCed.Scdsed.filt")
