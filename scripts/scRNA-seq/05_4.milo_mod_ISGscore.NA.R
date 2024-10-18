library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)


### 4.3 From Seurat object
library(Seurat)
set.seed(42)


## Load in the data
alldata_milo <- readRDS("mod/alldata_milo.rds")
alldata <- readRDS(file = "../harmony_seurat/Harmony_round1.rds")


# nhoods 1 に入っている細胞 nhoodsは2818個
str(alldata_milo@nhoods[,1])
##  num [1:66915] 0 0 0 0 0 0 0 0 0 0 ...

### id * nhoods
for (i in 1:2818) {
    if (i==1) {
        df <- cbind(rownames(colData(alldata_milo)), alldata_milo@nhoods[,i])
    } else {
        df <- cbind(df, alldata_milo@nhoods[,i])
    }
}

df <- as.data.frame(df)

# add colnames
colnames(df)[1] <- "Sample_Barcode"
for (i in 1:2818) {
    m=i+1
    colnames(df)[m] <- paste0("nhoods_", i)
}

# as.integer
for (i in 1:2818) {
    m=i+1
    df[,m] <- as.integer(df[,m])
}

### add columns case/ctrl
library(tibble)
casectrl <- alldata@meta.data %>% rownames_to_column(var = "Sample_Barcode") %>% select("Sample_Barcode", "SampleGroup")

df <- dplyr::left_join(df, casectrl, by="Sample_Barcode")

dir.create("ISG", recursive=T)
write.table(df, "ISG/SampleBarcodes_nhoods_casectrl.tsv", sep="\t", row.names=F, quote=F)



### calculate anti-viral ISG score in seurat
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
alldata <- AddModuleScore(
    object = alldata,
    features = antiviral_ISGs,
    ctrl = 100,
    name = 'antiviral_ISGs',
    seed=42
    )






df_isg <- data.frame(Nhood=NA, all_mean_ISG=NA, case_mean_ISG=NA, ctrl_mean_ISG=NA, Difference=NA, PValue=NA)

for (i in 1:2818) {
    m=i+1
    # case + ctrl
    all_ids <- df[df[,m]==1,]$Sample_Barcode
    # case
    case_ids <- df[df[,m]==1 & df$SampleGroup=="eHHV-6B",]$Sample_Barcode
    # ctrl
    ctrl_ids <- df[df[,m]==1 & df$SampleGroup=="Ctrl",]$Sample_Barcode
    
    # subset case+ctrl from seurat
    all_seurat <- subset(alldata, cells=all_ids)
    all_mean_isg <- mean(all_seurat@meta.data$antiviral_ISGs)
    # subset case from seurat
    if (length(case_ids)!=0) {
        case_seurat <- subset(alldata, cells=case_ids)
        case_mean_isg <- mean(case_seurat@meta.data$antiviral_ISGs)
    } else {
        case_mean_isg <- NA
    }
    # subset ctrl from seurat
    if (length(ctrl_ids)!=0) {
        ctrl_seurat <- subset(alldata, cells=ctrl_ids)
        ctrl_mean_isg <- mean(ctrl_seurat@meta.data$antiviral_ISGs)
    } else {
        ctrl_mean_isg <- NA
    }
    
    # Welch t test(対応のないパラメトリック検定)　←多重検定の観点から分散が均一でないことを確認しなくても第一選択という意見あり（等分散にも使える）
    if (length(case_ids)>1 & length(ctrl_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        welch_res <- t.test(case_seurat@meta.data$antiviral_ISGs, ctrl_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        ##          Welch Two Sample t-test
        ##  
        ##  data:  case_seurat@meta.data$antiviral_ISGs and ctrl_seurat@meta.data$antiviral_ISGs
        ##  t = 0.89671, df = 21.058, p-value = 0.38
        ##  alternative hypothesis: true difference in means is not equal to 0
        ##  95 percent confidence interval:
        ##   -0.04000583  0.10067732
        ##  sample estimates:
        ##   mean of x  mean of y
        ##  0.07878153 0.04844579
        PValue <- welch_res$p.value
        Difference = case_mean_isg - ctrl_mean_isg
    } else if ((length(case_ids)==1 & length(ctrl_ids)!=0) | (length(case_ids)!=0 & length(ctrl_ids)==1)) {
        Difference = case_mean_isg - ctrl_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        PValue <- NA
    } else {
        #PValue <- 1
        PValue <- NA
        Difference <- case_mean_isg - ctrl_mean_isg
    }
    df_isg[i,] <- c(i, all_mean_isg, case_mean_isg, ctrl_mean_isg, Difference, PValue)
}

da3_results_overlap5 <- read.csv("mod/da3_results_overlap5.csv", row.names=1)
Nhood_NhoodGroup <- da3_results_overlap5 %>% select('Nhood', 'NhoodGroup')

df_isg <- dplyr::left_join(df_isg, Nhood_NhoodGroup, by="Nhood")


write.table(df_isg, "ISG/df_isg_PvalueNA.tsv", sep="\t", quote=F, row.names=F)

# df_isg <- read.table("ISG/df_isg.tsv", sep="\t", header=T)
#df_isg[df_isg$PValue=="One",]$PValue <- ifelse(df_isg[df_isg$PValue=="One",]$Difference==0, 1, 0)
#df_isg[df_isg$PValue=="None",]$PValue <- 1
df_isg$PValue <- as.numeric(df_isg$PValue)

df_isg$SpatialFDR <- NA

## https://github.com/MarioniLab/miloR/issues/133
#df_isg$SpatialFDR <- miloR::graphSpatialFDR(
#                            x.nhoods=nhoods(alldata_milo),
#                            graph=miloR::graph(alldata_milo),
#                            weighting='k-distance',
#                            pvalues=df_isg$PValue,
#                            indices=nhoodIndex(alldata_milo),
#                            distances=nhoodDistances(alldata_milo),
#                            k=50)

df_isg[complete.cases(df_isg$PValue),]$SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$PValue),]$PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)


write.table(df_isg, "ISG/df_isg_PvalueNA.FDR.tsv", sep="\t", quote=F, row.names=F)


# from Seurat
UMAP_pl <- DimPlot(alldata, reduction = "umap", group.by = "predicted.celltype.l2") + NoLegend() +
             theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
UMAP_pl <- LabelClusters(plot = UMAP_pl, id = "predicted.celltype.l2", size=3)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(alldata_milo, df_isg, res_column="Difference", layout="UMAP", alpha=0.1)
p <- UMAP_pl + nh_graph_pl
ggsave(file="ISG/UMAP_Diff_alpha0.1.png", plot=p, width=10, height=5)

# FDR 1
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_isg, res_column="Difference", layout="UMAP", alpha=1)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG/UMAP_Diff_alpha1.png", plot=p, width=10, height=5)

# FDR 0.05
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_isg, res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG/UMAP_Diff_alpha0.05.png", plot=p, width=10, height=5)



### case_N or ctrl_N =0を除いてGroup内のcase平均値 vs ctrl平均値
# NAが含まれる行自体を無視して
# 正規性について　例えばGroup1のcase_mean_ISGでshapiro.testすると有意差ありで正規性なし　Group15でtestすると有意差なしで正規性あり
comp_df <- df_isg[complete.cases(df_isg),]
#comp_df <- df_isg[complete.cases(df_isg$Difference),]

# 対応のあるt検定
df_tres <- data.frame(NhoodGroup=NA, PValue=NA, CI_95=NA, Mean_of_Differences=NA)
for (i in  sort(comp_df$NhoodGroup) %>% unique()) {
    case_ISG <- comp_df[comp_df$NhoodGroup==i,]$case_mean_ISG
    ctrl_ISG <- comp_df[comp_df$NhoodGroup==i,]$ctrl_mean_ISG
    if (length(case_ISG) > 1 & length(ctrl_ISG) > 1) {
        t_res <- t.test(case_ISG, ctrl_ISG, paired=T)
        PValue=t_res$p.value
        CI_95=paste0(t_res$conf.int[1], "–", t_res$conf.int[2])
        Mean_of_Differences=t_res$estimate
    } else if (length(case_ISG)!=0 & length(ctrl_ISG)!=0) {
        PValue=0
        CI_95=NA
        Mean_of_Differences = mean(case_ISG) - mean(ctrl_ISG)
    } else {
        PValue=1
        CI_95=NA
        Mean_of_Differences=NA
    }
    df_tres[i,] <- c(i, PValue, CI_95, Mean_of_Differences)
}

df_tres$FDR <- p.adjust(df_tres$PValue, method="BH")
write.table(df_tres, "ISG/NhoodGroup_ttest.tsv", sep="\t", quote=F, row.names=F)

## r-ggplot2
df_isg <- read.table("ISG/df_isg.FDR.tsv", sep="\t", header=T)
comp_df <- df_isg[complete.cases(df_isg),]
df_plot <- reshape2::melt(data=comp_df, id.vars="NhoodGroup", measure.vars=c("ctrl_mean_ISG", "case_mean_ISG"))
colnames(df_plot) <- c("NhoodGroup", "variable", "mean_ISG")
library(dplyr)
df_plot <- df_plot %>% mutate(
    each_Nhood=case_when(
        variable=="case_mean_ISG" ~ "Case",
        variable=="ctrl_mean_ISG" ~ "Ctrl"))
df_plot$each_Nhood <- factor(df_plot$each_Nhood, levels=c("Ctrl", "Case"))
library(ggpubr)
library(ggsci)
p <- ggpaired(data=df_plot, x="each_Nhood", y="mean_ISG", color="each_Nhood") +
        stat_compare_means(comparisons=list(c("Ctrl", "Case")), label = "p.signif", method="t.test", paired=TRUE) +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup)
ggsave("ISG/NhoodGroup_meanISG.pdf", p, width=10, height=15)

### GroupSelection
df_plot_sub <- df_plot[df_plot$NhoodGroup==1 | df_plot$NhoodGroup==4 | df_plot$NhoodGroup==5 | df_plot$NhoodGroup==6 | df_plot$NhoodGroup==7 |
                        df_plot$NhoodGroup==8 | df_plot$NhoodGroup==9 | df_plot$NhoodGroup==10 | df_plot$NhoodGroup==12 | df_plot$NhoodGroup==13 |
                         df_plot$NhoodGroup==15 | df_plot$NhoodGroup==16 | df_plot$NhoodGroup==17 | df_plot$NhoodGroup==18 | df_plot$NhoodGroup==20,]
p <- ggpaired(data=df_plot_sub, x="each_Nhood", y="mean_ISG", color="each_Nhood") +
        stat_compare_means(comparisons=list(c("Ctrl", "Case")), label = "p.signif", method="t.test", paired=TRUE) +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup, ncol=8)
ggsave("ISG/NhoodGroup_meanISG_GroupSelection.pdf", p, width=10, height=7)


### levels=c("Case", "Ctrl")
df_plot$each_Nhood <- factor(df_plot$each_Nhood, levels=c("Case", "Ctrl"))
df_plot_sub <- df_plot[df_plot$NhoodGroup==1 | df_plot$NhoodGroup==4 | df_plot$NhoodGroup==5 | df_plot$NhoodGroup==6 | df_plot$NhoodGroup==7 |
                        df_plot$NhoodGroup==8 | df_plot$NhoodGroup==9 | df_plot$NhoodGroup==10 | df_plot$NhoodGroup==12 | df_plot$NhoodGroup==13 |
                         df_plot$NhoodGroup==15 | df_plot$NhoodGroup==16 | df_plot$NhoodGroup==17 | df_plot$NhoodGroup==18 | df_plot$NhoodGroup==20,]
p <- ggpaired(data=df_plot_sub, x="each_Nhood", y="mean_ISG", color="each_Nhood") +
        stat_compare_means(comparisons=list(c("Ctrl", "Case")), label = "p.signif", method="t.test", paired=TRUE) +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup, ncol=8)
ggsave("ISG/NhoodGroup_meanISG_GroupSelection.case-ctrl.pdf", p, width=10, height=7)
ggsave("ISG/NhoodGroup_meanISG_GroupSelection.case-ctrl.forPaper.pdf", p, width=10, height=5)

