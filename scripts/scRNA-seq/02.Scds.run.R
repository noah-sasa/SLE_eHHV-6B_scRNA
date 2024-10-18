library(dplyr)
library(tidyr)
library(tidyverse)
library(scds)
library(Seurat)
library(SingleCellExperiment)
set.seed(114)

args <- commandArgs(trailingOnly = T)
sampleid <- args[1]
dirname <- paste0("count_", sampleid)

group <- args[2]
institute <- args[3]


## Set up variables and parameters ##
out <- paste0("/mnt/e/scRNA_SLE/scds/", dirname, "/")
tenX_matrix <- paste0("/mnt/e/scRNA_SLE/cellranger/", dirname, "/outs/filtered_feature_bc_matrix/")

if (dir.exists(out)){
    unlink(out, recursive=T)
    dir.create(out)
} else {
    dir.create(out)
}

## Read in data
counts <- Read10X(as.character(tenX_matrix), gene.column = 2)

## Account for possibility that not just single cell data
if (is.list(counts)){
  sce <- SingleCellExperiment(list(counts=counts[[grep("Gene", names(counts))]]))
} else {
  sce <- SingleCellExperiment(list(counts=counts))
}

## Annotate doublet using binary classification based doublet scoring:
sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)

## Annotate doublet using co-expression based doublet scoring:
try({
    sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
})

### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
} else {
    print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
}

## Doublet scores are now available via colData:
colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType)
Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)

message("writing output")
write_delim(Doublets, paste0(out,"/scds_doublets_singlets.tsv"), "\t")


summary <- as.data.frame(table(Doublets$scds_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(out,"/scds_doublet_summary.tsv"), "\t")



# add sample and barcode to colData
colData(sce) <- colData(sce) %>%
    as.data.frame() %>%
    mutate(Sample=sampleid) %>%
    mutate(Barcode=rownames(colData(sce))) %>%
    mutate(SampleGroup=group) %>%
    mutate(Institute=institute) %>%
    DataFrame()


saveRDS(sce, file = paste0(out, "/sce"))

