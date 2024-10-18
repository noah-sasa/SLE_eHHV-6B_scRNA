### da3_results_overlap5
### DEG subset.row = NULL
dge <- testDiffExp(alldata_milo, da3_results_overlap5, design = ~ sex+scaled_age+group_id, meta.data = data.frame(colData(alldata_milo)),
                     subset.row = NULL, subset.nhoods=NULL)

for (i in 1:20) {
    dge_sub <- dge[[i]]
    dge_sub <- dge_sub[!is.na(dge_sub$adj.P.Val),]
    p <- paste0("p",i)
    library(ggthemes)
    dge_sub_thres <- dge_sub %>% 
                          mutate(threshold = adj.P.Val < 0.01 & abs(logFC) >= 0.58)
    dge_sub_thres$label <- NA
    dge_sub_thres$label[dge_sub_thres$threshold==TRUE] <- rownames(dge_sub_thres)[dge_sub_thres$threshold==TRUE]

    library(ggrepel)
    assign(p, ggplot(dge_sub_thres, aes(x=logFC,y=-log10(adj.P.Val), colour = threshold)) +
                    geom_point(alpha=0.5, size=2) +
                    ggtitle(paste0("Group",i)) +
                    xlab("log2 fold change") + 
                    ylab("-log10 adjusted p-value") +theme_foundation()+scale_colour_gdocs()   +
                    theme(legend.position = "none",
                        plot.background = element_blank(),
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) +
                    geom_text_repel(aes(label=label),
                        size=4,
                        box.padding = unit(0.35, "lines"),
                        point.pading = unit(0.3, "lines"),
                        max.overlaps = 20)
                    )

}
p <- p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20
ggsave(file="mod/DA3_overlap5_eacheGroup_DEG_ALL_Volcano.png", plot=p, width=20, height=20)


for (i in 1:20) {
    dge_sub <- dge[[i]]
    dge_sub <- dge_sub[!is.na(dge_sub$adj.P.Val),]
    p <- paste0("p",i)
    library(ggthemes)
    assign(p, ggplot(dge_sub, aes(x=P.Value)) +
                    ggtitle(paste0("Group",i)) +
                    geom_histogram(bins=50))

}
p <- p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20
ggsave(file="mod/DA3_overlap5_eacheGroup_DEG_ALL_PHistogram.png", plot=p, width=20, height=20)



for (i in 1:20) {
    dge_sub <- dge[[i]]
    write.csv(dge_sub, paste0("mod/DEG_100_5/DEG_ALL_Group",i,".csv"), quote=FALSE, row.names=TRUE)
}



### GroupのCtrl/eHHV_6Bの割合
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "group_id")



p <- ggplot(da3_results_overlap5, aes(x=group_id_fraction, fill=group_id)) + geom_histogram(bins=50) + facet_wrap(~ NhoodGroup) +
         scale_fill_gdocs() +
         theme(plot.background = element_blank())
ggsave(file="mod/DA3_overlap5_group_id_fraction.pdf", plot=p, width=7, height=7)

da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        eHHV_6B_fraction =
            if_else(
                group_id=="eHHV_6B",
                group_id_fraction,
                1 - group_id_fraction))


p <- ggplot(da3_results_overlap5, aes(x=eHHV_6B_fraction)) + geom_histogram(bins=50) + facet_wrap(~ NhoodGroup) +
         scale_fill_gdocs() +
         theme(plot.background = element_blank())
ggsave(file="mod/DA3_overlap5_eHHV_6B_fraction.pdf", plot=p, width=7, height=7)


### 積み上げ棒グラフ
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        Ctrl_fraction =
            if_else(
                group_id=="eHHV_6B",
                1 - group_id_fraction,
                group_id_fraction))

da3_results_overlap5$G_Nhood <- paste0(da3_results_overlap5$NhoodGroup, "_", da3_results_overlap5$Nhood)

library(dplyr)
df <- select(da3_results_overlap5, G_Nhood, eHHV_6B_fraction, Ctrl_fraction, NhoodGroup)
library(reshape2)
df <- melt(df)
head(df)
##    G_Nhood NhoodGroup         variable     value
##  1     1_1          1 eHHV_6B_fraction 0.3500000
##  2     2_2          2 eHHV_6B_fraction 0.4657534
##  3     3_3          3 eHHV_6B_fraction 0.1875000
##  4     4_4          4 eHHV_6B_fraction 0.3278689
##  5     4_5          4 eHHV_6B_fraction 0.1186441
##  6     5_6          5 eHHV_6B_fraction 0.2424242

library(plyr)
df <- ddply(df, .(G_Nhood),
                transform, pos = 100 - (cumsum(value*100) - (0.5 * value*100)))
head(df)
##    G_Nhood NhoodGroup         variable     value     pos
##  1     1_1          1 eHHV_6B_fraction 0.3500000 82.5000
##  2     1_1          1    Ctrl_fraction 0.6500000 32.5000
##  3  1_1014          1 eHHV_6B_fraction 0.2500000 87.5000
##  4  1_1014          1    Ctrl_fraction 0.7500000 37.5000
##  5  1_1028          1 eHHV_6B_fraction 0.4100719 79.4964
##  6  1_1028          1    Ctrl_fraction 0.5899281 29.4964

p <- ggplot() +
        geom_bar(data=df, aes(x=G_Nhood, y=value, fill=variable), stat="identity") + facet_wrap(~ NhoodGroup, scale="free") +
        #geom_text(data=df, aes(x=G_Nhood, y=pos, label=paste0(perc, "%")), size=4) +
        theme_classic() +
        #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.text.x = element_blank()) +
        scale_fill_gdocs() +
        scale_y_continuous(labels = function(x) paste0(x*100, "%"))

ggsave(plot=p, file="mod/DA3_overlap5_eHHV_6B_fraction.perc.pdf", height=7, width=7)





### sampleがG12に集中している可能性は？
### Groupのsample_id（Top1）の割合
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "group_id")

p <- ggplot(da3_results_overlap5, aes(x=sample_id_fraction, fill=sample_id)) + geom_histogram(bins=50) + facet_wrap(~ NhoodGroup) +
         theme(plot.background = element_blank())
ggsave(file="mod/DA3_overlap5_sample_id_fraction.pdf", plot=p, width=7, height=7)


### azimuth fraction
p <- ggplot(da3_results_overlap5, aes(x=predicted.celltype.l2_fraction, fill=predicted.celltype.l2)) + geom_histogram(bins=50) + facet_wrap(~ NhoodGroup) +
         theme(plot.background = element_blank())
ggsave(file="mod/DA3_overlap5_predicted.celltype.l2_fraction.pdf", plot=p, width=7, height=7)

p <- ggplot(da3_results_overlap5, aes(x=predicted.celltype.l2_fraction, fill=predicted.celltype.l2)) + geom_histogram(bins=50) + facet_wrap(~ NhoodGroup) + coord_cartesian(ylim=c(0,30)) +
         theme(plot.background = element_blank())
ggsave(file="mod/DA3_overlap5_predicted.celltype.l2_fraction.ylim30.pdf", plot=p, width=7, height=7)






#### 各サンプルの割合を無理やり出す
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE01")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE02")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE03")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE04")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE05")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE06")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE07")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE08")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE09")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE10")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE11")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE12")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE13")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE14")
da3_results_overlap5 <- annotateNhoods(alldata_milo, da3_results_overlap5, coldata_col = "SLE15")

da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE01_true_fraction =
            if_else(
                SLE01=="SLE01",
                SLE01_fraction,
                1 - SLE01_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE02_true_fraction =
            if_else(
                SLE02=="SLE02",
                SLE02_fraction,
                1 - SLE02_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE03_true_fraction =
            if_else(
                SLE03=="SLE03",
                SLE03_fraction,
                1 - SLE03_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE04_true_fraction =
            if_else(
                SLE04=="SLE04",
                SLE04_fraction,
                1 - SLE04_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE05_true_fraction =
            if_else(
                SLE05=="SLE05",
                SLE05_fraction,
                1 - SLE05_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE06_true_fraction =
            if_else(
                SLE06=="SLE06",
                SLE06_fraction,
                1 - SLE06_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE07_true_fraction =
            if_else(
                SLE07=="SLE07",
                SLE07_fraction,
                1 - SLE07_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE08_true_fraction =
            if_else(
                SLE08=="SLE08",
                SLE08_fraction,
                1 - SLE08_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE09_true_fraction =
            if_else(
                SLE09=="SLE09",
                SLE09_fraction,
                1 - SLE09_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE10_true_fraction =
            if_else(
                SLE10=="SLE10",
                SLE10_fraction,
                1 - SLE10_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE11_true_fraction =
            if_else(
                SLE11=="SLE11",
                SLE11_fraction,
                1 - SLE11_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE12_true_fraction =
            if_else(
                SLE12=="SLE12",
                SLE12_fraction,
                1 - SLE12_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE13_true_fraction =
            if_else(
                SLE13=="SLE13",
                SLE13_fraction,
                1 - SLE13_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE14_true_fraction =
            if_else(
                SLE14=="SLE14",
                SLE14_fraction,
                1 - SLE14_fraction))
da3_results_overlap5 <- da3_results_overlap5 %>%
    mutate(
        SLE15_true_fraction =
            if_else(
                SLE15=="SLE15",
                SLE15_fraction,
                1 - SLE15_fraction))


da3_results_overlap5$G_Nhood <- paste0(da3_results_overlap5$NhoodGroup, "_", da3_results_overlap5$Nhood)

library(dplyr)
df <- select(da3_results_overlap5, G_Nhood, SLE01_true_fraction, SLE02_true_fraction, SLE03_true_fraction, 
    SLE04_true_fraction, SLE05_true_fraction, SLE06_true_fraction, SLE07_true_fraction, SLE08_true_fraction, 
    SLE09_true_fraction, SLE10_true_fraction, SLE11_true_fraction, SLE12_true_fraction, SLE13_true_fraction, SLE14_true_fraction, 
    SLE15_true_fraction, NhoodGroup)
colnames(df) <- c("G_Nhood", "SLE01", "SLE02", "SLE03",
    "SLE04", "SLE05", "SLE06", "SLE07", "SLE08",
    "SLE09", "SLE10", "SLE11", "SLE12", "SLE13", "SLE14",
    "SLE15", "NhoodGroup")

library(reshape2)
df <- melt(df)

library(plyr)
df <- ddply(df, .(G_Nhood),
                transform, pos = 100 - (cumsum(value*100) - (0.5 * value*100)))

p <- ggplot() +
        geom_bar(data=df, aes(x=G_Nhood, y=value, fill=variable), stat="identity") + facet_wrap(~ NhoodGroup, scale="free") +
        #geom_text(data=df, aes(x=G_Nhood, y=pos, label=paste0(perc, "%")), size=4) +
        theme_classic() +
        #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.text.x = element_blank()) +
        #scale_fill_gdocs() +
        scale_y_continuous(labels = function(x) paste0(x*100, "%"))

ggsave(plot=p, file="mod/DA3_overlap5_sample_fraction.perc.pdf", height=7, width=7)



### in eHHV-6B
library(dplyr)
df <- select(da3_results_overlap5, G_Nhood, SLE01_true_fraction, 
    SLE02_true_fraction,
    SLE03_true_fraction,
    SLE04_true_fraction, NhoodGroup)

df$all_fraction <- df$SLE01_true_fraction + df$SLE02_true_fraction + df$SLE03_true_fraction + df$SLE04_true_fraction
df$SLE01_true02_fraction <- df$SLE01_true_fraction/df$all_fraction
df$SLE02_true02_fraction <- df$SLE02_true_fraction/df$all_fraction
df$SLE03_true02_fraction <- df$SLE03_true_fraction/df$all_fraction
df$SLE04_true02_fraction <- df$SLE04_true_fraction/df$all_fraction

df <- select(df, G_Nhood, SLE01_true02_fraction, 
    SLE02_true02_fraction,
    SLE03_true02_fraction,
    SLE04_true02_fraction, NhoodGroup)
colnames(df) <- c("G_Nhood", "SLE01",
    "SLE02",
    "SLE03",
    "SLE04", "NhoodGroup")
df$NhoodGroup <- as.factor(df$NhoodGroup)

df02 <- data.frame(df)

library(reshape2)
df <- melt(df)

library(plyr)
df <- ddply(df, .(G_Nhood),
                transform, pos = 100 - (cumsum(value*100) - (0.5 * value*100)))

p <- ggplot() +
        geom_bar(data=df, aes(x=G_Nhood, y=value, fill=variable), stat="identity") + facet_wrap(~ NhoodGroup, scale="free") +
        #geom_text(data=df, aes(x=G_Nhood, y=pos, label=paste0(perc, "%")), size=4) +
        theme_classic() +
        #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(axis.text.x = element_blank()) +
        scale_fill_gdocs() +
        scale_y_continuous(labels = function(x) paste0(x*100, "%"))

ggsave(plot=p, file="mod/DA3_overlap5_sample_eHHV-6B_fraction.perc.pdf", height=7, width=7)


### それぞれのGroupで平均Percentageが各症例がX%を超えず、かつ、Y%はあるようなGroupのみ解析したい
p <- ggplot() +
        geom_histogram(data=df, aes(x=value), fill="#b82e2e") + facet_grid(NhoodGroup ~ variable, scale="free") +
        #geom_text(data=df, aes(x=G_Nhood, y=pos, label=paste0(perc, "%")), size=4) +
                theme_foundation()+scale_colour_gdocs() +
            theme(plot.background = element_blank())
        #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        #scale_fill_gdocs() +
        #coord_cartesian(xlim=(c(0,1))) +
        #scale_x_continuous(labels = function(x) paste0(x*100, "%"))

ggsave(plot=p, file="mod/DA3_overlap5_sample_eHHV-6B_fraction.perc.histo.pdf", height=20, width=5)

df_mean <- data.frame(NhoodGroup=NA, SLE01=NA, SLE02=NA, SLE03=NA, SLE04=NA)
df_median <- data.frame(NhoodGroup=NA, SLE01=NA, SLE02=NA, SLE03=NA, SLE04=NA)
for (i in 1:20) {
    SLE01 <- mean(df[df$NhoodGroup==i & df$variable=="SLE01",]$value, na.rm=TRUE)
    SLE02 <- mean(df[df$NhoodGroup==i & df$variable=="SLE02",]$value, na.rm=TRUE)
    SLE03 <- mean(df[df$NhoodGroup==i & df$variable=="SLE03",]$value, na.rm=TRUE)
    SLE04 <- mean(df[df$NhoodGroup==i & df$variable=="SLE04",]$value, na.rm=TRUE)
    df_mean[i,] <- c(i, OI_00001, OI_00042, KR_00021, UO_00159)
    df_median[i,] <- c(i, OI_00001, OI_00042, KR_00021, UO_00159)
}


### 条件 最も高い中央値percentage >= 60% を除外
df_median[df_median$SLE01 < 0.60 & df_median$SLE02 < 0.60 & df_median$SLE03 < 0.60 & df_median$SLE04 < 0.60,]
##     NhoodGroup      SLE01      SLE02      SLE03      SLE04
##  1           1 0.33333333 0.34883721 0.22727273 0.04166667
##  4           4 0.21748311 0.18181818 0.49242424 0.00000000
##  5           5 0.21980676 0.54545455 0.10526316 0.00000000
##  6           6 0.35294118 0.07142857 0.00000000 0.50000000
##  7           7 0.20000000 0.17857143 0.00000000 0.57575758
##  8           8 0.05882353 0.46666667 0.37500000 0.00000000
##  9           9 0.18181818 0.19047619 0.55555556 0.00000000
##  10         10 0.33333333 0.28125000 0.26923077 0.03125000
##  12         12 0.55363985 0.23611111 0.16620690 0.00000000
##  13         13 0.41176471 0.35000000 0.14634146 0.05882353
##  15         15 0.25810811 0.41666667 0.03333333 0.00000000
##  16         16 0.58823529 0.04761905 0.20000000 0.14285714
##  17         17 0.58823529 0.15000000 0.21428571 0.00000000
##  18         18 0.20000000 0.23369565 0.01923077 0.50000000
##  20         20 0.00000000 0.33333333 0.33333333 0.33333333

# volcano again
dge <- testDiffExp(alldata_milo, da3_results_overlap5, design = ~ sex+scaled_age+group_id, meta.data = data.frame(colData(alldata_milo)),
                     subset.row = NULL, subset.nhoods=NULL)

for (i in 1:20) {
    dge_sub <- dge[[i]]
    dge_sub <- dge_sub[!is.na(dge_sub$adj.P.Val),]
    p <- paste0("p",i)
    library(ggthemes)
    dge_sub_thres <- dge_sub %>% 
                          mutate(threshold = adj.P.Val < 0.01 & abs(logFC) >= 0.58)
    dge_sub_thres$label <- NA
    dge_sub_thres$label[dge_sub_thres$threshold==TRUE] <- rownames(dge_sub_thres)[dge_sub_thres$threshold==TRUE]

    library(ggrepel)
    assign(p, ggplot(dge_sub_thres, aes(x=logFC,y=-log10(adj.P.Val), colour = threshold)) +
                    geom_point(alpha=0.5, size=2) +
                    ggtitle(paste0("Group",i)) +
                    xlab("log2 fold change") + 
                    ylab("-log10 adjusted p-value") +theme_foundation()+scale_colour_gdocs()   +
                    theme(legend.position = "none",
                        plot.background = element_blank(),
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) +
                    geom_text_repel(aes(label=label),
                        size=4,
                        box.padding = unit(0.35, "lines"),
                        point.pading = unit(0.3, "lines"),
                        max.overlaps = 20)
                    )

}
p <- p1+p4+p5+p6+p7+p8+p9+p10+p12+p13+p15+p16+p17+p18+p20+plot_layout(ncol=5)
ggsave(file="mod/DA3_overlap5_eacheGroup_DEG_ALL_Volcano_GroupSelection.png", plot=p, width=20, height=15)
ggsave(file="mod/DA3_overlap5_eacheGroup_DEG_ALL_Volcano_GroupSelection.pdf", plot=p, width=20, height=15)


### for Paper
### r-ggplot2
for (i in 1:20) {
    dge_sub <- read.table(paste0("mod/DEG_100_5/DEG_ALL_Group",i,".csv"), sep=",", header=TRUE)
    rownames(dge_sub) <- dge_sub$X
    dge_sub <- dge_sub[!is.na(dge_sub$adj.P.Val),]
    p <- paste0("p",i)
    library(ggthemes)
    library(tidyverse)
    dge_sub_thres <- dge_sub %>% 
                          mutate(threshold = adj.P.Val < 0.01 & abs(logFC) >= 0.58)
    dge_sub_thres$label <- NA
    dge_sub_thres$label[dge_sub_thres$threshold==TRUE] <- rownames(dge_sub_thres)[dge_sub_thres$threshold==TRUE]

    library(ggrepel)
    library(ggrastr)
    assign(p, ggplot(dge_sub_thres, aes(x=logFC,y=-log10(adj.P.Val), colour = threshold)) +
                    geom_point_rast(alpha=0.5, size=2) +
                    ggtitle(paste0("Group",i)) +
                    xlab("log2 fold change") + 
                    ylab("-log10 adjusted p-value") +theme_foundation()+scale_colour_gdocs()   +
                    theme(legend.position = "none",
                        plot.background = element_blank(),
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) +
                    geom_text_repel(aes(label=label),
                        size=4,
                        box.padding = unit(0.35, "lines"),
                        point.pading = unit(0.3, "lines"),
                        max.overlaps = 20)
                    )
}
library(patchwork)
p <- p1+p4+p5+p6+p7+p8+p9+p10+p12+p13+p15+p16+p17+p18+p20+plot_layout(ncol=5)
ggsave(file="mod/DA3_overlap5_eacheGroup_DEG_ALL_Volcano_GroupSelection.forPaper.pdf", plot=p, width=20, height=15)





write.csv(da3_results_overlap5, "mod/da3_results_overlap5.csv", quote=FALSE, row.names=TRUE)
