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

### immunosuppressants
alldata@meta.data$immunosuppressants <- NA
alldata@meta.data[alldata@meta.data$Sample=="SLE01",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE02",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE03",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE04",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE05",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE06",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE07",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE08",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE09",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE10",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE11",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE12",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE13",]$immunosuppressants <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE14",]$immunosuppressants <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE15",]$immunosuppressants <- "Yes"

### prednisolone
alldata@meta.data$prednisolone <- NA
alldata@meta.data[alldata@meta.data$Sample=="SLE01",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE02",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE03",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE04",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE05",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE06",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE07",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE08",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE09",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE10",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE11",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE12",]$prednisolone <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE13",]$prednisolone <- "No"
alldata@meta.data[alldata@meta.data$Sample=="SLE14",]$prednisolone <- "Yes"
alldata@meta.data[alldata@meta.data$Sample=="SLE15",]$prednisolone <- "Yes"


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
casectrl <- alldata@meta.data %>% rownames_to_column(var = "Sample_Barcode") %>% select("Sample_Barcode", "SampleGroup", "immunosuppressants", "prednisolone")

df <- dplyr::left_join(df, casectrl, by="Sample_Barcode")

dir.create("ISG_immnosuppressants", recursive=T)
write.table(df, "ISG_immnosuppressants/SampleBarcodes_nhoods_casectrl.tsv", sep="\t", row.names=F, quote=F)



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






df_isg <- data.frame(Nhood=NA, all_mean_ISG=NA, all_immuno_mean_ISG=NA, all_NOimmuno_mean_ISG=NA, all_Difference=NA, all_PValue=NA,
    case_immuno_mean_ISG=NA, case_NOimmuno_mean_ISG=NA, ctrl_immuno_mean_ISG=NA, ctrl_NOimmuno_mean_ISG=NA,
    case_immuno_case_NOimmuno_Difference=NA, case_immuno_ctrl_immuno_Difference=NA, case_immuno_ctrl_NOimmuno_Difference=NA,
    case_NOimmuno_ctrl_immuno_Difference=NA, case_NOimmuno_ctrl_NOimmuno_Difference=NA, ctrl_immuno_ctrl_NOimmuno_Difference=NA,
    case_immuno_case_NOimmuno_PValue=NA, case_immuno_ctrl_immuno_PValue=NA, case_immuno_ctrl_NOimmuno_PValue=NA,
    case_NOimmuno_ctrl_immuno_PValue=NA, case_NOimmuno_ctrl_NOimmuno_PValue=NA, ctrl_immuno_ctrl_NOimmuno_PValue=NA)

for (i in 1:2818) {
    m=i+1
    # case + ctrl
    all_ids <- df[df[,m]==1,]$Sample_Barcode
    # immuno
    immuno_ids <- df[df[,m]==1 & df$immunosuppressants=="Yes",]$Sample_Barcode
    # NOimmuno
    NOimmuno_ids <- df[df[,m]==1 & df$immunosuppressants=="No",]$Sample_Barcode
    # case_immuno
    case_immuno_ids <- df[df[,m]==1 & df$immunosuppressants=="Yes" & df$SampleGroup=="eHHV-6B",]$Sample_Barcode
    # case_NOimmuno
    case_NOimmuno_ids <- df[df[,m]==1 & df$immunosuppressants=="No" & df$SampleGroup=="eHHV-6B",]$Sample_Barcode
    # ctrl_immuno
    ctrl_immuno_ids <- df[df[,m]==1 & df$immunosuppressants=="Yes" & df$SampleGroup=="Ctrl",]$Sample_Barcode
    # ctrl_NOimmuno
    ctrl_NOimmuno_ids <- df[df[,m]==1 & df$immunosuppressants=="No" & df$SampleGroup=="Ctrl",]$Sample_Barcode
    
    # subset case+ctrl from seurat
    all_seurat <- subset(alldata, cells=all_ids)
    all_mean_isg <- mean(all_seurat@meta.data$antiviral_ISGs)
    # subset immuno from seurat
    if (length(immuno_ids)!=0) {
        immuno_seurat <- subset(alldata, cells=immuno_ids)
        immuno_mean_isg <- mean(immuno_seurat@meta.data$antiviral_ISGs)
    } else {
        immuno_mean_isg <- NA
    }
    # subset NOimmuno from seurat
    if (length(NOimmuno_ids)!=0) {
        NOimmuno_seurat <- subset(alldata, cells=NOimmuno_ids)
        NOimmuno_mean_isg <- mean(NOimmuno_seurat@meta.data$antiviral_ISGs)
    } else {
        NOimmuno_mean_isg <- NA
    }
    # subset case_immuno from seurat
    if (length(case_immuno_ids)!=0) {
        case_immuno_seurat <- subset(alldata, cells=case_immuno_ids)
        case_immuno_mean_isg <- mean(case_immuno_seurat@meta.data$antiviral_ISGs)
    } else {
        case_immuno_mean_isg <- NA
    }
    # subset case_NOimmuno from seurat
    if (length(case_NOimmuno_ids)!=0) {
        case_NOimmuno_seurat <- subset(alldata, cells=case_NOimmuno_ids)
        case_NOimmuno_mean_isg <- mean(case_NOimmuno_seurat@meta.data$antiviral_ISGs)
    } else {
        case_NOimmuno_mean_isg <- NA
    }
    # subset ctrl_immuno from seurat
    if (length(ctrl_immuno_ids)!=0) {
        ctrl_immuno_seurat <- subset(alldata, cells=ctrl_immuno_ids)
        ctrl_immuno_mean_isg <- mean(ctrl_immuno_seurat@meta.data$antiviral_ISGs)
    } else {
        ctrl_immuno_mean_isg <- NA
    }
    # subset ctrl_NOimmuno from seurat
    if (length(ctrl_NOimmuno_ids)!=0) {
        ctrl_NOimmuno_seurat <- subset(alldata, cells=ctrl_NOimmuno_ids)
        ctrl_NOimmuno_mean_isg <- mean(ctrl_NOimmuno_seurat@meta.data$antiviral_ISGs)
    } else {
        ctrl_NOimmuno_mean_isg <- NA
    }
    
    # Welch t test(対応のないパラメトリック検定)　←多重検定の観点から分散が均一でないことを確認しなくても第一選択という意見あり（等分散にも使える）
    if (length(immuno_ids)>1 & length(NOimmuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        all_welch_res <- t.test(immuno_seurat@meta.data$antiviral_ISGs, NOimmuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
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
        all_PValue <- all_welch_res$p.value
        all_Difference = immuno_mean_isg - NOimmuno_mean_isg
    } else if ((length(immuno_ids)==1 & length(NOimmuno_ids)!=0) | (length(immuno_ids)!=0 & length(NOimmuno_ids)==1)) {
        all_Difference = immuno_mean_isg - NOimmuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        all_PValue <- NA
    } else {
        #PValue <- 1
        all_PValue <- NA
        all_Difference <- immuno_mean_isg - NOimmuno_mean_isg
    }

    if (length(case_immuno_ids)>1 & length(case_NOimmuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        case_immuno_case_NOimmuno_welch_res <- t.test(case_immuno_seurat@meta.data$antiviral_ISGs, case_NOimmuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        case_immuno_case_NOimmuno_PValue <- case_immuno_case_NOimmuno_welch_res$p.value
        case_immuno_case_NOimmuno_Difference = case_immuno_mean_isg - case_NOimmuno_mean_isg
    } else if ((length(case_immuno_ids)==1 & length(case_NOimmuno_ids)!=0) | (length(case_immuno_ids)!=0 & length(case_NOimmuno_ids)==1)) {
        case_immuno_case_NOimmuno_Difference = case_immuno_mean_isg - case_NOimmuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        case_immuno_case_NOimmuno_PValue <- NA
    } else {
        #PValue <- 1
        case_immuno_case_NOimmuno_PValue <- NA
        case_immuno_case_NOimmuno_Difference <- case_immuno_mean_isg - case_NOimmuno_mean_isg
    }

    if (length(case_immuno_ids)>1 & length(ctrl_immuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        case_immuno_ctrl_immuno_welch_res <- t.test(case_immuno_seurat@meta.data$antiviral_ISGs, ctrl_immuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        case_immuno_ctrl_immuno_PValue <- case_immuno_ctrl_immuno_welch_res$p.value
        case_immuno_ctrl_immuno_Difference = case_immuno_mean_isg - ctrl_immuno_mean_isg
    } else if ((length(case_immuno_ids)==1 & length(ctrl_immuno_ids)!=0) | (length(case_immuno_ids)!=0 & length(ctrl_immuno_ids)==1)) {
        case_immuno_ctrl_immuno_Difference = case_immuno_mean_isg - ctrl_immuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        case_immuno_ctrl_immuno_PValue <- NA
    } else {
        #PValue <- 1
        case_immuno_ctrl_immuno_PValue <- NA
        case_immuno_ctrl_immuno_Difference <- case_immuno_mean_isg - ctrl_immuno_mean_isg
    }

    if (length(case_immuno_ids)>1 & length(ctrl_NOimmuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        case_immuno_ctrl_NOimmuno_welch_res <- t.test(case_immuno_seurat@meta.data$antiviral_ISGs, ctrl_NOimmuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        case_immuno_ctrl_NOimmuno_PValue <- case_immuno_ctrl_NOimmuno_welch_res$p.value
        case_immuno_ctrl_NOimmuno_Difference = case_immuno_mean_isg - ctrl_NOimmuno_mean_isg
    } else if ((length(case_immuno_ids)==1 & length(ctrl_NOimmuno_ids)!=0) | (length(case_immuno_ids)!=0 & length(ctrl_NOimmuno_ids)==1)) {
        case_immuno_ctrl_NOimmuno_Difference = case_immuno_mean_isg - ctrl_NOimmuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        case_immuno_ctrl_NOimmuno_PValue <- NA
    } else {
        #PValue <- 1
        case_immuno_ctrl_NOimmuno_PValue <- NA
        case_immuno_ctrl_NOimmuno_Difference <- case_immuno_mean_isg - ctrl_NOimmuno_mean_isg
    }

    if (length(case_NOimmuno_ids)>1 & length(ctrl_immuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        case_NOimmuno_ctrl_immuno_welch_res <- t.test(case_NOimmuno_seurat@meta.data$antiviral_ISGs, ctrl_immuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        case_NOimmuno_ctrl_immuno_PValue <- case_NOimmuno_ctrl_immuno_welch_res$p.value
        case_NOimmuno_ctrl_immuno_Difference = case_NOimmuno_mean_isg - ctrl_immuno_mean_isg
    } else if ((length(case_NOimmuno_ids)==1 & length(ctrl_immuno_ids)!=0) | (length(case_NOimmuno_ids)!=0 & length(ctrl_immuno_ids)==1)) {
        case_NOimmuno_ctrl_immuno_Difference = case_NOimmuno_mean_isg - ctrl_immuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        case_NOimmuno_ctrl_immuno_PValue <- NA
    } else {
        #PValue <- 1
        case_NOimmuno_ctrl_immuno_PValue <- NA
        case_NOimmuno_ctrl_immuno_Difference <- case_NOimmuno_mean_isg - ctrl_immuno_mean_isg
    }

    if (length(case_NOimmuno_ids)>1 & length(ctrl_NOimmuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        case_NOimmuno_ctrl_NOimmuno_welch_res <- t.test(case_NOimmuno_seurat@meta.data$antiviral_ISGs, ctrl_NOimmuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        case_NOimmuno_ctrl_NOimmuno_PValue <- case_NOimmuno_ctrl_NOimmuno_welch_res$p.value
        case_NOimmuno_ctrl_NOimmuno_Difference = case_NOimmuno_mean_isg - ctrl_NOimmuno_mean_isg
    } else if ((length(case_NOimmuno_ids)==1 & length(ctrl_NOimmuno_ids)!=0) | (length(case_NOimmuno_ids)!=0 & length(ctrl_NOimmuno_ids)==1)) {
        case_NOimmuno_ctrl_NOimmuno_Difference = case_NOimmuno_mean_isg - ctrl_NOimmuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        case_NOimmuno_ctrl_NOimmuno_PValue <- NA
    } else {
        #PValue <- 1
        case_NOimmuno_ctrl_NOimmuno_PValue <- NA
        case_NOimmuno_ctrl_NOimmuno_Difference <- case_NOimmuno_mean_isg - ctrl_NOimmuno_mean_isg
    }

    if (length(ctrl_immuno_ids)>1 & length(ctrl_NOimmuno_ids)>1) {
        # t.testは2:2以上　1:多の場合は平均値が大きいか小さいか、一緒かは自明
        ctrl_immuno_ctrl_NOimmuno_welch_res <- t.test(ctrl_immuno_seurat@meta.data$antiviral_ISGs, ctrl_NOimmuno_seurat@meta.data$antiviral_ISGs, alternative="two.sided", paired=FALSE, var.equal=FALSE, conf.level=0.95)
        ctrl_immuno_ctrl_NOimmuno_PValue <- ctrl_immuno_ctrl_NOimmuno_welch_res$p.value
        ctrl_immuno_ctrl_NOimmuno_Difference = ctrl_immuno_mean_isg - ctrl_NOimmuno_mean_isg
    } else if ((length(ctrl_immuno_ids)==1 & length(ctrl_NOimmuno_ids)!=0) | (length(ctrl_immuno_ids)!=0 & length(ctrl_NOimmuno_ids)==1)) {
        ctrl_immuno_ctrl_NOimmuno_Difference = ctrl_immuno_mean_isg - ctrl_NOimmuno_mean_isg
        #PValue <- ifelse(Difference==0, 1, 0)
        ctrl_immuno_ctrl_NOimmuno_PValue <- NA
    } else {
        #PValue <- 1
        ctrl_immuno_ctrl_NOimmuno_PValue <- NA
        ctrl_immuno_ctrl_NOimmuno_Difference <- ctrl_immuno_mean_isg - ctrl_NOimmuno_mean_isg
    }


    df_isg[i,] <- c(i, all_mean_isg, immuno_mean_isg, NOimmuno_mean_isg, all_Difference, all_PValue,
        case_immuno_mean_isg, case_NOimmuno_mean_isg, ctrl_immuno_mean_isg, ctrl_NOimmuno_mean_isg,
        case_immuno_case_NOimmuno_Difference, case_immuno_ctrl_immuno_Difference, case_immuno_ctrl_NOimmuno_Difference,
        case_NOimmuno_ctrl_immuno_Difference, case_NOimmuno_ctrl_NOimmuno_Difference, ctrl_immuno_ctrl_NOimmuno_Difference,
        case_immuno_case_NOimmuno_PValue, case_immuno_ctrl_immuno_PValue, case_immuno_ctrl_NOimmuno_PValue,
        case_NOimmuno_ctrl_immuno_PValue, case_NOimmuno_ctrl_NOimmuno_PValue, ctrl_immuno_ctrl_NOimmuno_PValue)
}

da3_results_overlap5 <- read.csv("mod/da3_results_overlap5.csv", row.names=1)
Nhood_NhoodGroup <- da3_results_overlap5 %>% select('Nhood', 'NhoodGroup')

df_isg <- dplyr::left_join(df_isg, Nhood_NhoodGroup, by="Nhood")


write.table(df_isg, "ISG_immnosuppressants/df_isg_PvalueNA.tsv", sep="\t", quote=F, row.names=F)

# df_isg <- read.table("ISG/df_isg.tsv", sep="\t", header=T)
#df_isg[df_isg$PValue=="One",]$PValue <- ifelse(df_isg[df_isg$PValue=="One",]$Difference==0, 1, 0)
#df_isg[df_isg$PValue=="None",]$PValue <- 1
df_isg$all_PValue <- as.numeric(df_isg$all_PValue)
df_isg$case_immuno_case_NOimmuno_PValue <- as.numeric(df_isg$case_immuno_case_NOimmuno_PValue)
df_isg$case_immuno_ctrl_immuno_PValue <- as.numeric(df_isg$case_immuno_ctrl_immuno_PValue)
df_isg$case_immuno_ctrl_NOimmuno_PValue <- as.numeric(df_isg$case_immuno_ctrl_NOimmuno_PValue)
df_isg$case_NOimmuno_ctrl_immuno_PValue <- as.numeric(df_isg$case_NOimmuno_ctrl_immuno_PValue)
df_isg$case_NOimmuno_ctrl_NOimmuno_PValue <- as.numeric(df_isg$case_NOimmuno_ctrl_NOimmuno_PValue)
df_isg$ctrl_immuno_ctrl_NOimmuno_PValue <- as.numeric(df_isg$ctrl_immuno_ctrl_NOimmuno_PValue)

df_isg$all_SpatialFDR <- NA
df_isg$case_immuno_case_NOimmuno_SpatialFDR <- NA
df_isg$case_immuno_ctrl_immuno_SpatialFDR <- NA
df_isg$case_immuno_ctrl_NOimmuno_SpatialFDR <- NA
df_isg$case_NOimmuno_ctrl_immuno_SpatialFDR <- NA
df_isg$case_NOimmuno_ctrl_NOimmuno_SpatialFDR <- NA
df_isg$ctrl_immuno_ctrl_NOimmuno_SpatialFDR <- NA

df_isg[complete.cases(df_isg$all_PValue),]$all_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$all_PValue),]$all_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$case_immuno_case_NOimmuno_PValue),]$case_immuno_case_NOimmuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$case_immuno_case_NOimmuno_PValue),]$case_immuno_case_NOimmuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$case_immuno_ctrl_immuno_PValue),]$case_immuno_ctrl_immuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$case_immuno_ctrl_immuno_PValue),]$case_immuno_ctrl_immuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$case_immuno_ctrl_NOimmuno_PValue),]$case_immuno_ctrl_NOimmuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$case_immuno_ctrl_NOimmuno_PValue),]$case_immuno_ctrl_NOimmuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$case_NOimmuno_ctrl_immuno_PValue),]$case_NOimmuno_ctrl_immuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$case_NOimmuno_ctrl_immuno_PValue),]$case_NOimmuno_ctrl_immuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$case_NOimmuno_ctrl_NOimmuno_PValue),]$case_NOimmuno_ctrl_NOimmuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$case_NOimmuno_ctrl_NOimmuno_PValue),]$case_NOimmuno_ctrl_NOimmuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)
df_isg[complete.cases(df_isg$ctrl_immuno_ctrl_NOimmuno_PValue),]$ctrl_immuno_ctrl_NOimmuno_SpatialFDR <- miloR::graphSpatialFDR(
                            x.nhoods=nhoods(alldata_milo),
                            graph=miloR::graph(alldata_milo),
                            weighting='k-distance',
                            pvalues=df_isg[complete.cases(df_isg$ctrl_immuno_ctrl_NOimmuno_PValue),]$ctrl_immuno_ctrl_NOimmuno_PValue,
                            indices=nhoodIndex(alldata_milo),
                            distances=nhoodDistances(alldata_milo),
                            k=50)


write.table(df_isg, "ISG_immnosuppressants/df_isg_PvalueNA.FDR.tsv", sep="\t", quote=F, row.names=F)


# from Seurat
UMAP_pl <- DimPlot(alldata, reduction = "umap", group.by = "predicted.celltype.l2") + NoLegend() +
             theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
UMAP_pl <- LabelClusters(plot = UMAP_pl, id = "predicted.celltype.l2", size=3)

## Plot neighbourhood graph
# FDR 0.05
library(tidyverse)
df_tmp <- df_isg %>% rename("SpatialFDR" = all_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = all_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_all_Diff_alpha0.05.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = case_immuno_case_NOimmuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = case_immuno_case_NOimmuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_case_NOimmuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_case_NOimmuno_Diff_alphaNA.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = case_immuno_ctrl_immuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = case_immuno_ctrl_immuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_ctrl_immuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_ctrl_immuno_Diff_alpha0.05.pdf", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_ctrl_immuno_Diff_alphaNA.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = case_immuno_ctrl_NOimmuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = case_immuno_ctrl_NOimmuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_ctrl_NOimmuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_immuno_ctrl_NOimmuno_Diff_alphaNA.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = case_NOimmuno_ctrl_immuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = case_NOimmuno_ctrl_immuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_NOimmuno_ctrl_immuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_NOimmuno_ctrl_immuno_Diff_alphaNA.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = case_NOimmuno_ctrl_NOimmuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = case_NOimmuno_ctrl_NOimmuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_NOimmuno_ctrl_NOimmuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_case_NOimmuno_ctrl_NOimmuno_Diff_alphaNA.png", plot=p, width=10, height=5)

df_tmp <- df_isg %>% rename("SpatialFDR" = ctrl_immuno_ctrl_NOimmuno_SpatialFDR)
df_tmp <- df_tmp %>% rename("Difference" = ctrl_immuno_ctrl_NOimmuno_Difference)
nh_graph_pl1 <- plotNhoodGraphDA(alldata_milo, df_tmp[complete.cases(df_tmp$SpatialFDR),], res_column="Difference", layout="UMAP", alpha=0.05)
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_ctrl_immuno_ctrl_NOimmuno_Diff_alpha0.05.png", plot=p, width=10, height=5)
alldata_milo_tmp <- alldata_milo
colData(alldata_milo_tmp)["Difference"] <- NA
colData(alldata_milo_tmp)[unlist(nhoodIndex(alldata_milo_tmp)[df_tmp$Nhood]),"Difference"] <- df_tmp[,"Difference"]
nh_graph_pl1 <- plotNhoodGraph(alldata_milo_tmp, colour_by = "Difference", layout="UMAP")
p <- UMAP_pl + nh_graph_pl1
ggsave(file="ISG_immnosuppressants/UMAP_ctrl_immuno_ctrl_NOimmuno_Diff_alphaNA.png", plot=p, width=10, height=5)



df_plot <- reshape2::melt(data=df_isg, id.vars="Nhood", measure.vars=c("all_immuno_mean_ISG", "all_NOimmuno_mean_ISG", "case_immuno_mean_ISG", "case_NOimmuno_mean_ISG", "ctrl_immuno_mean_ISG", "ctrl_NOimmuno_mean_ISG"))
df_plot <- dplyr::left_join(df_plot, Nhood_NhoodGroup, by="Nhood")
colnames(df_plot) <- c("Nhood", "variable", "mean_ISG", "NhoodGroup")

library(ggplot2)
library(ggsci)
p <- ggplot(data=df_plot, aes(x=variable, y=mean_ISG, color=variable)) +
        geom_boxplot() +
        geom_point() +
        geom_line() +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup)
ggsave("ISG_immnosuppressants/NhoodGroup_meanISG.pdf", p, width=10, height=15)


### case_immuno vs ctrl_immuno
df_isg_immuno <- df_isg %>% select('Nhood', 'NhoodGroup', 'case_immuno_mean_ISG', 'ctrl_immuno_mean_ISG')
comp_df <- df_isg_immuno[complete.cases(df_isg_immuno),]
# 対応のあるt検定
df_tres <- data.frame(NhoodGroup=NA, PValue=NA, CI_95=NA, Mean_of_Differences=NA)
for (i in  sort(comp_df$NhoodGroup) %>% unique()) {
    case_ISG <- comp_df[comp_df$NhoodGroup==i,]$case_immuno_mean_ISG
    ctrl_ISG <- comp_df[comp_df$NhoodGroup==i,]$ctrl_immuno_mean_ISG
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

df_tres$BH <- p.adjust(df_tres$PValue, method="BH")
write.table(df_tres, "ISG_immnosuppressants/case_immuno_ctrl_immuno.NhoodGroup_ttest.tsv", sep="\t", quote=F, row.names=F)

## r-ggplot2
df_isg <- read.table("ISG_immnosuppressants/df_isg_PvalueNA.FDR.tsv", sep="\t", header=T)
df_isg_immuno <- df_isg %>% select('Nhood', 'NhoodGroup', 'case_immuno_mean_ISG', 'ctrl_immuno_mean_ISG')
comp_df <- df_isg_immuno[complete.cases(df_isg_immuno),]
df_plot <- reshape2::melt(data=comp_df, id.vars="NhoodGroup", measure.vars=c("case_immuno_mean_ISG", "ctrl_immuno_mean_ISG"))
colnames(df_plot) <- c("NhoodGroup", "variable", "mean_ISG")
library(dplyr)
df_plot <- df_plot %>% mutate(
    each_Nhood=case_when(
        variable=="case_immuno_mean_ISG" ~ "Case",
        variable=="ctrl_immuno_mean_ISG" ~ "Ctrl"))
df_plot$each_Nhood <- factor(df_plot$each_Nhood, levels=c("Case", "Ctrl"))
library(ggpubr)
library(ggsci)
p <- ggpaired(data=df_plot, x="each_Nhood", y="mean_ISG", color="each_Nhood") +
        stat_compare_means(comparisons=list(c("Ctrl", "Case")), label = "p.signif", method="t.test", paired=TRUE) +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup)
ggsave("ISG_immnosuppressants/case_immuno_ctrl_immuno.NhoodGroup_meanISG.pdf", p, width=10, height=15)

### GroupSelection
df_plot_sub <- df_plot[df_plot$NhoodGroup==1 | df_plot$NhoodGroup==4 | df_plot$NhoodGroup==5 | df_plot$NhoodGroup==6 | df_plot$NhoodGroup==7 |
                        df_plot$NhoodGroup==8 | df_plot$NhoodGroup==9 | df_plot$NhoodGroup==10 | df_plot$NhoodGroup==12 | df_plot$NhoodGroup==13 |
                         df_plot$NhoodGroup==15 | df_plot$NhoodGroup==16 | df_plot$NhoodGroup==17 | df_plot$NhoodGroup==18 | df_plot$NhoodGroup==20,]
p <- ggpaired(data=df_plot_sub, x="each_Nhood", y="mean_ISG", color="each_Nhood") +
        stat_compare_means(comparisons=list(c("Ctrl", "Case")), label = "p.signif", method="t.test", paired=TRUE) +
        scale_color_lancet() +
        facet_wrap(~NhoodGroup, ncol=8)
ggsave("ISG_immnosuppressants/case_immuno_ctrl_immuno.NhoodGroup_meanISG_GroupSelection.case-ctrl.pdf", p, width=10, height=7)





