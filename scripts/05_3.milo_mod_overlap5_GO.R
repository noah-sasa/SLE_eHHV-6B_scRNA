library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(ggplot2)



### lFC 1.5 p < 0.01
get_ReactomePA <- function(i){
    library(ReactomePA)
    df <- read.table(paste0("../../milo/mod/DEG_100_5/DEG_ALL_Group",i,".csv"), sep=",", header=T)
    eg <- clusterProfiler::bitr(df$X, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    colnames(eg) <- c("X", "ENTREZID")
    
    df <- dplyr::inner_join(df, eg, by="X")

    geneList <- df[df$adj.P.Val < 0.01, ]$logFC                     # adj.P.Val < 0.01
    names(geneList) <- df[df$adj.P.Val < 0.01, ]$ENTREZID            # adj.P.Val < 0.01
    geneList = sort(geneList, decreasing = TRUE)

    de <- names(geneList)[abs(geneList) > 0.58]     # logFC < -1.5 | logFC > 1.5

    x <- enrichPathway(gene=de, pvalueCutoff=0.05, readable=T)
    dir.create("enrichPathway", recursive=TRUE)
    saveRDS(x, paste0("enrichPathway/enrichPathway_Group", i, ".lFC.adjP0.011.5.rds"))
    write.csv(as.data.frame(x), paste0("enrichPathway/enrichPathway_Group", i, ".lFC1.5.adjP0.01.csv"), quote=FALSE, row.names=TRUE)

    library(ggplot2)
    library(ggrepel)

    pbar <- barplot(x, showCategory=8) + theme_classic()

    pdot <- dotplot(x, showCategory=15) + theme_classic()

    x2 <- enrichplot::pairwise_termsim(x)
    pemap <- emapplot(x2)

    pcnet <- cnetplot(x2, categorySize='pvalue', foldChange=geneList)


    y <- gsePathway(geneList, nPerm=10000,
                            pvalueCutoff=0.2,
                            pAdjustMethod="BH", verbose=FALSE)

    res <- as.data.frame(y)
    dir.create("gsePathway", recursive=TRUE)
    saveRDS(y, paste0("gsePathway/gsePathway_Group", i, ".lFC1.5.adjP0.01.rds"))
    write.csv(res, paste0("gsePathway/gsePathway_Group", i, ".lFC1.5.adjP0.01.csv"), quote=FALSE, row.names=TRUE)

    y2 <- enrichplot::pairwise_termsim(y)
    pemap_gsea <- emapplot(y2, color="pvalue")

    #pgsea <- gseaplot(y2, geneSetID = head(y$ID, n=1)) + theme_classic()

    #pviewPathway <- viewPathway(head(y$Description, n=1), readable=TRUE, foldChange=geneList) + theme_classic()

    library(patchwork)

    p <- pbar + pdot + pemap + pcnet + pemap_gsea + plot_layout(ncol = 5)
    dir.create("figures", recursive=TRUE)
    ggsave(paste0("figures/Group_", i, "_matome.lFC1.5.adjP0.01.pdf"), p, width=50, height=10, limitsize = FALSE)
    pmatome <- paste0("pmatome_", i)
    assign(pmatome, pbar + pdot + pemap + pcnet + pemap_gsea)
}

for (i in 1:20) {
    tryCatch({
        get_ReactomePA(i)
    }, error=function(e){})
}


de_list <- list()
gene_list <- list()
for (i in 1:20) {
    library(ReactomePA)
    df <- read.table(paste0("../../milo/mod/DEG_100_5/DEG_ALL_Group",i,".csv"), sep=",", header=T)
    eg <- clusterProfiler::bitr(df$X, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    colnames(eg) <- c("X", "ENTREZID")
    
    df <- dplyr::inner_join(df, eg, by="X")
    gene_list[[i]] <- df$ENTREZID

    geneList <- df[df$adj.P.Val < 0.01, ]$logFC                     # adj.P.Val < 0.01
    names(geneList) <- df[df$adj.P.Val < 0.01, ]$ENTREZID            # adj.P.Val < 0.01
    geneList = sort(geneList, decreasing = TRUE)

    de <- names(geneList)[abs(geneList) > 0.58]     # logFC < -1.5 | logFC > 1.5

    de_list[[i]] <- de
}


### GroupSelection 1 4 5 6 7 8 9 10 12 13 15 16 17 18 20
de_list_sub <- de_list[c(1,4,5,6,7,8,9,10,12,13,15,16,17,18,20)]
gene_list_sub <- gene_list[c(1,4,5,6,7,8,9,10,12,13,15,16,17,18,20)]
# GO CC
res <- compareCluster(de_list_sub, fun="enrichGO",
                    universe=gene_list_sub,
                            OrgDb=org.Hs.eg.db,
                            ont="CC",               # cellular component
                    pAdjustMethod="BH",
                    pvalueCutoff=0.01,
                    qvalueCutoff=0.05,
                            )
p <- dotplot(res) + theme_classic()
ggsave("Comparing_enriched_GO_CC.lFC1.5.adjP0.01_GroupSelection.pdf", p, width=5, height=5)
write.csv(res, "Comparing_enriched_GO_CC.lFC1.5.adjP0.01_GroupSelection.csv", quote=FALSE, row.names=FALSE)

# GO BP
res <- compareCluster(de_list_sub, fun="enrichGO",
                    universe=gene_list_sub,
                            OrgDb=org.Hs.eg.db,
                            ont="BP",               # biological process
                    pAdjustMethod="BH",
                    pvalueCutoff=0.01,
                    qvalueCutoff=0.05,
                            )
p <- dotplot(res) + theme_classic()
ggsave("Comparing_enriched_GO_BP.lFC1.5.adjP0.01_GroupSelection.pdf", p, width=7, height=7)
write.csv(res, "Comparing_enriched_GO_BP.lFC1.5.adjP0.01_GroupSelection.csv", quote=FALSE, row.names=FALSE)

# GO MF
res <- compareCluster(de_list_sub, fun="enrichGO",
                    universe=gene_list_sub,
                            OrgDb=org.Hs.eg.db,
                            ont="MF",               # molecular function
                    pAdjustMethod="BH",
                    pvalueCutoff=0.01,
                    qvalueCutoff=0.05,
                            )
p <- dotplot(res) + theme_classic()
ggsave("Comparing_enriched_GO_MF.lFC1.5.adjP0.01_GroupSelection.pdf", p, width=5, height=5)
write.csv(res, "Comparing_enriched_GO_MF.lFC1.5.adjP0.01_GroupSelection.csv", quote=FALSE, row.names=FALSE)

# GO ALL
res <- compareCluster(de_list_sub, fun="enrichGO",
                    universe=gene_list_sub,
                            OrgDb=org.Hs.eg.db,
                            ont="all",               # all: CC+BP+MF
                    pAdjustMethod="BH",
                    pvalueCutoff=0.01,
                    qvalueCutoff=0.05,
                            )
p <- dotplot(res, split="ONTOLOGY") + theme_classic() + facet_grid(ONTOLOGY~., scale="free")
ggsave("Comparing_enriched_GO_ALL.lFC1.5.adjP0.01_GroupSelection.pdf", p, width=7, height=14)
write.csv(res, "Comparing_enriched_GO_ALL.lFC1.5.adjP0.01_GroupSelection.csv", quote=FALSE, row.names=FALSE)





### GO BP
get_GO_BP <- function(i){
    library(clusterProfiler)
    df <- read.table(paste0("../../milo/mod/DEG_100_5/DEG_ALL_Group",i,".csv"), sep=",", header=T)
    eg <- clusterProfiler::bitr(df$X, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    colnames(eg) <- c("X", "ENTREZID")
    
    df <- dplyr::inner_join(df, eg, by="X")
    allgene <- df$ENTREZID

    geneList <- df[df$adj.P.Val < 0.01, ]$logFC                     # adj.P.Val < 0.01
    names(geneList) <- df[df$adj.P.Val < 0.01, ]$ENTREZID            # adj.P.Val < 0.01
    geneList = sort(geneList, decreasing = TRUE)

    de <- names(geneList)[abs(geneList) > 0.58]     # FC < -1.5 | FC > 1.5

    x <- enrichGO(gene=de,
                    universe=allgene,
                    OrgDb=org.Hs.eg.db,
                    ont="BP",
                    pAdjustMethod="BH",
                    pvalueCutoff=0.01,
                    qvalueCutoff=0.05,
                    readable=T)
    dir.create("enrichGO", recursive=TRUE)
    saveRDS(x, paste0("enrichGO/enrichGO_BP_Group", i, ".lFC1.5.adjP0.01.rds"))
    write.csv(as.data.frame(x), paste0("enrichGO/enrichGO_BP_Group", i, ".lFC1.5.adjP0.01.csv"), quote=FALSE, row.names=TRUE)

    library(ggplot2)
    library(ggrepel)

    pbar <- barplot(x, showCategory=12) + theme_classic()

    pdot <- dotplot(x, showCategory=15) + theme_classic()

    x2 <- enrichplot::pairwise_termsim(x)
    pemap <- emapplot(x2)

    pcnet <- cnetplot(x2, categorySize='pvalue', foldChange=geneList)

    pgo <- goplot(x2)

    pheat <- heatplot(x2, foldChange=geneList)


    y <- gseGO(geneList, nPerm=10000,
                            OrgDb=org.Hs.eg.db,
                            ont="BP",
                            minGSSize=120,
                            pvalueCutoff=0.05,
                            pAdjustMethod="BH", verbose=FALSE)

    res <- as.data.frame(y)
    dir.create("gseGO", recursive=TRUE)
    saveRDS(y, paste0("gseGO/gseGO_BP_Group", i, ".lFC1.5.adjP0.01.rds"))
    write.csv(res, paste0("gseGO/gseGO_BP_Group", i, ".lFC1.5.adjP0.01.csv"), quote=FALSE, row.names=TRUE)

    y2 <- enrichplot::pairwise_termsim(y)
    pemap_gsea <- emapplot(y2, color="pvalue")

    #pridge <- ridgeplot(y2, showCategory=12)       # need "ggridges"

    #pgsea <- gseaplot(y2, geneSetID = head(y$ID, n=1)) + theme_classic()

    #pviewPathway <- viewPathway(head(y$Description, n=1), readable=TRUE, foldChange=geneList) + theme_classic()

    library(patchwork)

    p <- pbar + pdot + pemap + pcnet + pgo + pheat + pemap_gsea + plot_layout(ncol = 7)
    dir.create("figuresGO", recursive=TRUE)
    ggsave(paste0("figuresGO/GO_BP_Group_", i, "_matome.lFC1.5.adjP0.01.pdf"), p, width=70, height=10, limitsize = FALSE)
    #ggsave(paste0("figuresGO/GP_BP_Group_", i, "_Heatmap.lFC1.5.adjP0.01.pdf"), pheat, width=10, height=10, limitsize = FALSE)
    pmatome <- paste0("pmatome_", i)
    assign(pmatome, pbar + pdot + pemap + pcnet + pgo + pemap_gsea)
}

for (i in 1:20) {
    tryCatch({
        get_GO_BP(i)
    }, error=function(e){})
}

# i=18 GSEAで0 enriched terms found
# p <- pbar + pdot + pemap + pcnet + pgo + pheat + plot_layout(ncol = 6)
# ggsave(paste0("figuresGO/GO_BP_Group_", i, "_matome.lFC1.5.adjP0.01.pdf"), p, width=60, height=10, limitsize = FALSE)

