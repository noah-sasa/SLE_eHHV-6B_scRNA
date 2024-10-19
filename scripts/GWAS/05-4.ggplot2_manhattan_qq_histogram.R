### https://sites.google.com/view/stuck-in-the-shallow-end/home/generate-manhattan-plots-with-ggplot2-and-ggrastr
### https://slowkow.com/notes/ggplot2-qqplot/


library(data.table)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(tidyverse)
library(ggrepel)

args=commandArgs(trailingOnly = T)
maf = args[1]
covariate = args[2]
ADDorDOM = args[3]
tsv=paste0("data03_maf", maf, "/plink2_logistic_results_", covariate, ".PHENO1.glm.logistic.", ADDorDOM, "only.removeNA.tsv")

df <- as.data.frame(fread(tsv, head=TRUE, sep="\t"), stringsAsFactors=FALSE)
df$CHROM <- factor(df$CHROM, levels=as.character(1:22))

## Calculate cummulative value of position. This will be used as the position on X-axis on Manhattan plot.
df$POS2 = df$POS
for (CHROM in 2:22) {
    tmp = CHROM - 1
    tmp2 = max(as.numeric(df[df$CHROM==tmp, "POS2"]), na.rm=TRUE)
    df[df$CHROM==CHROM, "POS2"] = df[df$CHROM==CHROM, "POS"] + tmp2
}

## Calculate center positions of each chromosome on X-axis. This will be used to set the positions of chromosome numbers.
df2 <- data.frame(CHROM=1:22, center=NA, stringsAsFactors=FALSE)
for (CHROM in 1:22) {
    start = min(df[df$CHROM==CHROM, "POS2"], na.rm=TRUE)
    end = max(df[df$CHROM==CHROM, "POS2"], na.rm=TRUE)
    center = mean(c(start, end))
    df2[df2$CHROM==CHROM, "center"] <- center
}

### Make a repeating pattern of two colors. My favorite is:
colour1="blue4"
colour2="orange3"

### Get Bonferroni-corrected threshold.
bonf=0.05/nrow(df)
genomewide=5e-8
suggestive=1e-5

### Highlight top SNPs
annotatePval <- 0.01
topHits = subset(df, P <= annotatePval)

# extract top snps in each chromosome
topHits <- topHits[order(topHits$P),]
topSNPs <- NULL

for (i in unique(topHits$CHROM)) {
    chrSNPs <- topHits[topHits$CHR == i,]
    topSNPs <- rbind(topSNPs, chrSNPs[1,])
}

topSNPs$label <- topSNPs$ID
label <- topSNPs %>% select(ID, label)

df <- left_join(df, label, by='ID')



### Manhattan plot:
g1 <- ggplot(df, aes(x=POS2, y=-log10(P), colour=CHROM, label=label)) +
        geom_point_rast(size=0.5) +
        #geom_hline(yintercept=-log10(bonf), linetype="dashed", color="red") +
        geom_hline(yintercept=-log10(genomewide), color="red") +
        geom_hline(yintercept=-log10(suggestive), color="blue") +
        geom_text_repel(size=1.5, color='black',
                        force_pull   = 0, 
                        nudge_x = 0.5,
                        box.padding = 0.5,
                        nudge_y = 0.5,
                        min.segment.length = 0, # draw all lines no matter how short
                        segment.size = 0.2,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 45,
                        ) +
        scale_colour_manual(values=rep(c(colour1, colour2), times=11)) +
        theme_classic() +
        theme(legend.position="NONE") +
        scale_x_continuous(breaks=df2$center, labels=c(1:16,"",18,"",20,"",22)) +
        xlab("Chromosome") +
        ylab(expression(paste(-log[10],"(",italic(P),")")))
ggsave(plot=g1, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".Manhattan.pdf"), width=10, height=4)

### for Paper
gg1 <- ggplot(df, aes(x=POS2, y=-log10(P), colour=CHROM, label=label)) +
        geom_point_rast(size=0.5) +
        geom_point(data=df[df$P<genomewide,], aes(x=POS2, y=-log10(P), colour=CHROM, label=label), size=1.5) +
        #geom_hline(yintercept=-log10(bonf), linetype="dashed", color="red") +
        geom_hline(yintercept=-log10(genomewide), color="red") +
        geom_hline(yintercept=-log10(suggestive), color="blue") +
        geom_text_repel(size=1.5, color='black',
                        force_pull   = 0, 
                        nudge_x = 0.5,
                        box.padding = 0.5,
                        nudge_y = 0.5,
                        min.segment.length = 0, # draw all lines no matter how short
                        segment.size = 0.2,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 45,
                        ) +
        scale_colour_manual(values=rep(c(colour1, colour2), times=11)) +
        theme_classic() +
        theme(legend.position="NONE") +
        scale_x_continuous(breaks=df2$center, labels=c(1:16,"",18,"",20,"",22)) +
        xlab("Chromosome") +
        ylab(expression(paste(-log[10],"(",italic(P),")")))
ggsave(plot=gg1, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".Manhattan.forPaper.pdf"), width=10, height=2.5)



### Lambda: a measure of inflated p-value
chisq <- qchisq(1-df$P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)

### QQ plot:
df_nrow = nrow(df)
exp.pval = (1:df_nrow -0.5)/df_nrow
exp.pval.log = as.data.frame(-log10(exp.pval))
var.pval = df$P[order(df$P)]
var.pval.log = as.data.frame(-log10(var.pval))
N = df_nrow
cupper = -log10(qbeta(0.95, 1:N, N-1:N+1))
clower = -log10(qbeta(1-0.95, 1:N, N-1:N+1))
df3 = cbind(exp.pval.log, var.pval.log, cupper, clower)
colnames(df3) = c("expected", "var", "cup", "clow")
g2 <- ggplot(df3) +
        geom_abline(slope=1, intercept=0, color='grey') +
        #geom_line(aes(expected, cup), linetype=2) +
        #geom_line(aes(expected, clow), linetype=2) +
        geom_point_rast(aes(x=expected, y=var), color='black', size=0.5) +
        theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank()) +
        xlab(expression(paste(-log[10], "(expected ", italic(P), ")"))) +
        ylab(expression(paste(-log[10], "(observed ", italic(P), ")"))) +
        theme_classic() +
        annotate(geom="text",
            x=-Inf,
            y=Inf,
            hjust=-0.15,
            vjust=1+0.15*3,
            label=sprintf("lambda = %.2f", lambda),
            size = 6)
ggsave(plot=g2, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".QQ.pdf"), width=4, height=4)

g2a <- ggplot(df3) +
        geom_abline(slope=1, intercept=0, color='grey') +
        geom_line(aes(expected, cup), linetype=2) +
        geom_line(aes(expected, clow), linetype=2) +
        geom_point_rast(aes(x=expected, y=var), color='black', size=0.5) +
        theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank()) +
        xlab(expression(paste(-log[10], "(expected ", italic(P), ")"))) +
        ylab(expression(paste(-log[10], "(observed ", italic(P), ")"))) +
        theme_classic() +
        annotate(geom="text",
            x=-Inf,
            y=Inf,
            hjust=-0.15,
            vjust=1+0.15*3,
            label=sprintf("lambda = %.2f", lambda),
            size = 6)
ggsave(plot=g2a, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".QQ_ci.pdf"), width=4, height=4)

### Histogram of p-values
g3 <- ggplot(df) +
        geom_histogram(aes(x=P), bins=25, color='white', size=0.3, boundary=0.5) +
        theme_classic() +
        xlab("") +
        ylab("") +
        ggtitle(expression(paste("Histogram of ", italic(P), "-values"))) +
        scale_y_continuous(expand = c(0.02, 0)) +
        theme(axis.line.x = element_line(size=0.5),
                axis.ticks.x = element_line(size=0.5),
                panel.grid.major.y = element_line(size=0.5),
                panel.grid.minor.y = element_line(size=0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
ggsave(plot=g3, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".P_histogram.pdf"), width=6, height=4)



### merge
layout <- "
AAABB
CCCCC
"
g = g3 + g2 + g1 + plot_layout(design = layout) + plot_annotation(tag_levels='a') + theme(plot.tag = element_text(size=24))
ggsave(plot=g, file=paste0("data03_maf", maf, "/ggplot2_logistic.", covariate, ".", ADDorDOM, ".pdf"), width=10, height=8)


