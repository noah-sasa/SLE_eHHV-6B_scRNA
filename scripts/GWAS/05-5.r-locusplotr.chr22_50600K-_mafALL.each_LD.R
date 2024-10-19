# conda-pack locusplotr

args=commandArgs(trailingOnly = T)
maf = args[1]
covariate = args[2]
ADDorDOM = args[3]



library(tidyverse)
library(data.table)
# lead snps rs1253523660, rs1278832979, rs56207913, rs76038028, rs146190586
for (lead_rsid in c("rs1253523660", "rs1278832979", "rs56207913", "rs76038028", "rs146190586")) {
    ld <- paste0(lead_rsid, "_ld")
    assign(ld, fread(paste0("/work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/chr22_T2T/", lead_rsid, "_LDscore.ld"), header=TRUE))
}
rs1253523660_ld <- rs1253523660_ld %>% dplyr::select(SNP_B, R2)
colnames(rs1253523660_ld) <- c("SNP_B", "rs1253523660")
rs1278832979_ld <- rs1278832979_ld %>% dplyr::select(SNP_B, R2)
colnames(rs1278832979_ld) <- c("SNP_B", "rs1278832979")
rs56207913_ld <- rs56207913_ld %>% dplyr::select(SNP_B, R2)
colnames(rs56207913_ld) <- c("SNP_B", "rs56207913")
rs76038028_ld <- rs76038028_ld %>% dplyr::select(SNP_B, R2)
colnames(rs76038028_ld) <- c("SNP_B", "rs76038028")
rs146190586_ld <- rs146190586_ld %>% dplyr::select(SNP_B, R2)
colnames(rs146190586_ld) <- c("SNP_B", "rs146190586")

ld <- full_join(rs1253523660_ld, rs1278832979_ld, by="SNP_B")
ld <- full_join(ld, rs56207913_ld, by="SNP_B")
ld <- full_join(ld, rs76038028_ld, by="SNP_B")
ld <- full_join(ld, rs146190586_ld, by="SNP_B")

ld_melt <- reshape2::melt(ld)
colnames(ld_melt) <- c("SNP_B", "lead_rsid", "R2")
ld_melt <- ld_melt[!is.na(ld_melt$R2),]
max_ld <- ld_melt %>% group_by(SNP_B) %>% filter(R2==max(R2))   # https://stackoverflow.com/questions/43720750/getting-the-name-of-the-variable-with-maximum-value-in-a-molten-data-frame

# 同じR2を取るlead_snpがあると重複してしまう
max_ld %>% group_by(SNP_B) %>% filter(n()>1)
#   # A tibble: 8 × 3
#   # Groups:   SNP_B [4]
#     SNP_B                          lead_rsid         R2
#     <chr>                          <fct>          <dbl>
#   1 chr22_50830883_CTTTTTTT_CTTTTT rs1253523660 0.00618
#   2 chr22_51322107_G_T             rs1253523660 0.0933
#   3 chr22_51322166_GAGGGTT_TAGGGTT rs1253523660 0.0349
#   4 chr22_51322798_TTGGGG_*        rs1253523660 0.00204
#   5 chr22_50830883_CTTTTTTT_CTTTTT rs146190586  0.00618
#   6 chr22_51322107_G_T             rs146190586  0.0933
#   7 chr22_51322166_GAGGGTT_TAGGGTT rs146190586  0.0349
#   8 chr22_51322798_TTGGGG_*        rs146190586  0.00204

#max_ld <- max_ld[-which((max_ld$SNP_B=="chr22_50830883_CTTTTTTT_CTTTTT" | max_ld$SNP_B=="chr22_51322107_G_T" | max_ld$SNP_B=="chr22_51322166_GAGGGTT_TAGGGTT" | max_ld$SNP_B=="chr22_51322798_TTGGGG_*") & max_ld$lead_rsid=="rs1253523660"),]

# 末端lead_rsid優先で
max_ld <- max_ld %>% distinct(SNP_B, .keep_all=TRUE)






df <- fread(paste0("/work23/home/nsasa/data/CHM13v2.0_hhv6/gwas04_sampleploidy/data03_maf0/plink2_logistic_results_", covariate, ".PHENO1.glm.logistic.", ADDorDOM, "only.removeNA.tsv"), sep="\t", header=TRUE)

df
##           CHROM      POS                                   ID A1
##        1:     1     3897 chr1_3897_CGCCTGCTGGCAGCTGAGGACACT_C  C
##        2:     1     5183               chr1_5183_TATCTCTTAG_T  T
##        3:     1     5245                       chr1_5245_CT_C  C
##        4:     1    22395                         rs1308908416 TA
##        5:     1    22831                    chr1_22831_CCTA_C  C
##       ---
##  5796999:    22 51281850                          rs371535564  T
##  5797000:    22 51287415                          rs147159805  A
##  5797001:    22 51287944                            rs2238837  C
##  5797002:    22 51294075                           rs28729663  A
##  5797003:    22 51300942                          rs200589453  T
##                                 A2       OR LOG(OR)_SE    Z_STAT           P
##        1: CGCCTGCTGGCAGCTGAGGACACT 1.226950   0.627795  0.325792 0.744582000
##        2:               TATCTCTTAG 0.154762   0.754717 -2.472280 0.013425600
##        3:                       CT 2.042330   0.598627  1.192880 0.232916000
##        4:                        T 1.672190   0.562254  0.914420 0.360496000
##        5:                     CCTA 0.303777   1.042810 -1.142540 0.253228000
##       ---
##  5796999:                       TA 1.378390   0.423665  0.757485 0.448759000
##  5797000:                        G 0.192944   1.439750 -1.142810 0.253119000
##  5797001:                        A 2.677760   0.341618  2.883280 0.003935580
##  5797002:                        G 0.261649   0.724418 -1.850800 0.064198600
##  5797003:                     TAAG 4.168010   0.379583  3.760550 0.000169543
##           ERRCODE
##        1:       .
##        2:       .
##        3:       .
##        4:       .
##        5:       .
##       ---
##  5796999:       .
##  5797000:       .
##  5797001:       .
##  5797002:       .
##  5797003:       .

#lead_rsid="rs76038028"
df[df$ID==lead_rsid,]
##     CHROM      POS         ID A1 A2    OR LOG(OR)_SE  Z_STAT           P ERRCODE
##  1:    22 51132381 rs76038028  A  T 15.92      0.505 5.48035 4.24486e-08       .
lead_chrom=df[df$ID==lead_rsid,]$CHROM
lead_pos=df[df$ID==lead_rsid,]$POS

# extract 1Mb
#   plot_distance = 500000
#   plot_start = lead_pos - plot_distance
#   plot_end = lead_pos + plot_distance

plot_start = 50600000
plot_end = 51400000

df_locus <- df[df$CHROM==lead_chrom & (df$POS>=plot_start & df$POS<=plot_end),]

df_locus
##        CHROM      POS                  ID A1   A2       OR LOG(OR)_SE    Z_STAT
##     1:    22 50632644           rs9617086  A    T 1.098760   0.365902  0.257405
##     2:    22 50632680            rs138240  T    G 1.085270   0.366841  0.223052
##     3:    22 50632737           rs9617033  A    G 1.098760   0.365902  0.257405
##     4:    22 50632879           rs9617087  T    C 1.139730   0.363097  0.360218
##     5:    22 50632882 chr22_50632882_A_AT AT    A 1.226770   0.352607  0.579629
##    ---
##  1700:    22 51281850         rs371535564  T   TA 1.378390   0.423665  0.757485
##  1701:    22 51287415         rs147159805  A    G 0.192944   1.439750 -1.142810
##  1702:    22 51287944           rs2238837  C    A 2.677760   0.341618  2.883280
##  1703:    22 51294075          rs28729663  A    G 0.261649   0.724418 -1.850800
##  1704:    22 51300942         rs200589453  T TAAG 4.168010   0.379583  3.760550
##                  P ERRCODE
##     1: 0.796866000       .
##     2: 0.823495000       .
##     3: 0.796866000       .
##     4: 0.718684000       .
##     5: 0.562165000       .
##    ---
##  1700: 0.448759000       .
##  1701: 0.253119000       .
##  1702: 0.003935580       .
##  1703: 0.064198600       .
##  1704: 0.000169543       .


#   ### LD data about rs76038028 from T2T_phased 1KG&Pangenome
#   ld <- fread(paste0("/work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/chr22_T2T/", lead_rsid, "_LDscore.ld"), header=TRUE)
#   #ld <- fread("/work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/chr22_phased_T2T/rs76038028_LDscore.ld", header=TRUE)
#   
#   ld
#   ##         CHR_A     BP_A      SNP_A CHR_B     BP_B
#   ##      1:    22 51132381 rs76038028    22 50132537
#   ##      2:    22 51132381 rs76038028    22 50132615
#   ##      3:    22 51132381 rs76038028    22 50132657
#   ##      4:    22 51132381 rs76038028    22 50132815
#   ##      5:    22 51132381 rs76038028    22 50132854
#   ##     ---
#   ##  22262:    22 51132381 rs76038028    22 51323242
#   ##  22263:    22 51132381 rs76038028    22 51323248
#   ##  22264:    22 51132381 rs76038028    22 51323254
#   ##  22265:    22 51132381 rs76038028    22 51323259
#   ##  22266:    22 51132381 rs76038028    22 51323263
#   ##                                                      SNP_B          R2        DP
#   ##      1:                                        rs557967053 5.40152e-05 1.0000000
#   ##      2:                                         rs35123758 8.79593e-05 0.0112363
#   ##      3:                                        rs560775536 5.40152e-05 1.0000000
#   ##      4:                                       rs1925932039 5.40152e-05 1.0000000
#   ##      5:                                        rs562445549 5.40152e-05 1.0000000
#   ##     ---
#   ##  22262:         chr22_51323242_G_GGGGTTGGGGTTGGGGTTGGGGTTA 5.40152e-05 1.0000000
#   ##  22263: chr22_51323248_GGGGTTGGGGTTGGGGT_AGGGTTGGGGTTGGGGT 1.08580e-04 1.0000000
#   ##  22264:                                 chr22_51323254_G_A 1.62368e-04 1.0000000
#   ##  22265:                                chr22_51323259_TG_T 5.44585e-05 1.0000000
#   ##  22266:                                chr22_51323263_GT_G 4.28571e-01 1.0000000
#   
#   # left_join
#   df_locus_r2 <- left_join(df_locus, ld, by=c('ID' = 'SNP_B'))

df_locus_r2 <- left_join(df_locus, max_ld, by=c('ID' = 'SNP_B'))


### recombination_rate
recomb <- fread("/work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/recomb1000GAvg/recomb1000GAvg.chr22.hg38_to_chm13v2.wig", header=FALSE) 
colnames(recomb) <- c("CHROM", "start", "end", "recomb_rate")

recomb
##         CHROM     start       end recomb_rate
##      1:    22   5890193  11473636    0.000000
##      2:    22  17063184  17063215    0.000000
##      3:    22   6337594  15870986    0.000000
##      4:    22 132167288 132207908    0.000000
##      5:    22  78030963  78031150    0.000000
##     ---
##  43432:    22  51287944  51291633    0.759015
##  43433:    22  51291633  51292203    0.759396
##  43434:    22  51292203  51297169    0.765501
##  43435:    22  51297169  51298706    0.776718
##  43436:    22  51298706  51304902    0.000000

# genetic_map
##         chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
##      1:  22 15870986            0.0000000     0.000000000
##      2:  22 15891057            0.0000000     0.000000000
##      3:  22 15903792            0.0000000     0.000000000
##      4:  22 15908158            0.0000000     0.000000000
##      5:  22 15919300            0.8914971     0.009893835
##     ---
##  43279:  22 51291633            0.7590147    73.922965294
##  43280:  22 51292203            0.7593958    73.923398150
##  43281:  22 51297169            0.7655011    73.927199629
##  43282:  22 51298706            0.7767175    73.928393443
##  43283:  22 51304902            0.0000000    73.928393443

# extract 1Mb
recomb_zoom <- recomb[recomb$start>=plot_start & recomb$end<=plot_end,]

recomb_zoom
##       CHROM    start      end recomb_rate
##    1:    22 50632644 50632680   0.0158627
##    2:    22 50632680 50633034   0.0685531
##    3:    22 50633034 50633249   1.6270600
##    4:    22 50633249 50633813   0.5624200
##    5:    22 50633813 50633866   0.5624200
##   ---
##  461:    22 51287944 51291633   0.7590150
##  462:    22 51291633 51292203   0.7593960
##  463:    22 51292203 51297169   0.7655010
##  464:    22 51297169 51298706   0.7767180
##  465:    22 51298706 51304902   0.0000000




# Create color codes and labels
df_locus_r2 <- df_locus_r2 %>%
    #mutate(color_code = as.character(
    #    cut(
    #        as.numeric(R2),
    #        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #        labels = c("blue4", "skyblue", "darkgreen", "orange", "red"),
    #        include.lowest = TRUE)
    #    )) %>%
    mutate(color_code = case_when(
        lead_rsid=="rs1253523660" ~ "#034788",
        lead_rsid=="rs1278832979" ~ "#E60808",
        lead_rsid=="rs56207913" ~ "#44B543",
        lead_rsid=="rs76038028" ~ "#0698B1",
        lead_rsid=="rs146190586" ~ "#915F9E")
    ) %>%
    #mutate(legend_label = as.character(
    #    cut(
    #        as.numeric(R2),
    #        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #        labels = c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1"),
    #        include.lowest = TRUE)
    #    )) %>%
    mutate(lead = case_when(
        ID == lead_rsid ~ TRUE,
        TRUE ~ FALSE)
    ) %>%
    mutate(label = case_when(
        ID == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
        )) %>%
    #mutate(color_code = case_when(
    #    ID == lead_rsid ~ "purple",
    #    TRUE ~ color_code
    #    )) %>%
    #mutate(color_code = forcats::fct_expand(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
    #mutate(color_code = forcats::fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
    mutate(color_code = forcats::fct_expand(color_code, "#034788", "#E60808", "#44B543", "#0698B1", "#915F9E", NA)) %>%
    mutate(color_code = forcats::fct_relevel(color_code, "#034788", "#E60808", "#44B543", "#0698B1", "#915F9E", NA)) %>%
    #mutate(legend_label = case_when(
    #    ID == lead_rsid ~ "Ref",
    #    TRUE ~ legend_label
    #    )) %>%
    #mutate(legend_label = forcats::fct_expand(legend_label, "Ref", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2")) %>%
    #mutate(legend_label = forcats::fct_relevel(legend_label, "Ref", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"))
    mutate(R2 = case_when(
        is.na(R2) ~ -log10(P)/10,
        TRUE ~ R2)
    )


# Make plot (sample non-significant p-values to reduce overplotting)
plot_pvalue_threshold = 0.1

# https://notchained.hatenablog.com/entry/2016/10/02/011204
scale_to_value1 <- function(values) scales::rescale(values, to = range(-log10(df_locus_r2$P)))
#scale_to_value2 <- function(values) scales::rescale(values, to = range(recomb_zoom$recomb_rate))                                       # 2nd y axisの最大値をrecomb_rateの最大値に
scale_to_value2 <- function(values) scales::rescale(values, to = range(recomb_zoom$recomb_rate / (max(recomb_zoom$recomb_rate)/100)))   # 2nd y axisの最大値を100に(locuszoom準拠)

library(ggplot2)
library(ggrastr)
regional_assoc_plot <- df_locus_r2 %>%
            distinct(ID, .keep_all = TRUE) %>%
            filter(P < plot_pvalue_threshold | R2 > 0.2) %>% # improve pverplotting
            bind_rows(df_locus_r2 %>%
                        filter(P >= plot_pvalue_threshold & R2 < 0.2)) %>%
            arrange(desc(color_code)) %>%
            ggplot() +
            geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
            #geom_step(data=recomb_zoom, aes(x=end, y=scale_to_value1(recomb_rate)), color='blue', na.rm=TRUE) +                            # 2nd y axisの最大値をrecomb_rateの最大値に
            geom_step(data=recomb_zoom, aes(x=end, y=scale_to_value1(recomb_rate) * (max(recomb_rate)/100)), color='blue', na.rm=TRUE) +    # 2nd y axisの最大値を100に(locuszoom準拠)
            geom_point_rast(aes(POS, -log10(P), fill = factor(color_code), size = lead, alpha = R2, shape = lead), color='black', stroke=0.5) +
            ggrepel::geom_label_repel(aes(POS, -log10(P), label = label),
                                        size = 4,
                                        color = "black",
                                        fontface = "bold",
                                        fill = "white",
                                        min.segment.length = 0,
                                        box.padding = 1,
                                        alpha = 1,
                                        nudge_y = 4
                                        ) +
            scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = c(levels(forcats::fct_drop(df_locus_r2$lead_rsid)), 'No data'), na.translate = FALSE) +
            scale_alpha_identity(parse(text = "r^2"), guide = "legend") +
            scale_size_manual(values = c(3, 5), guide = "none") +
            scale_shape_manual(values = c(21, 23), guide = "none") +
            scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            scale_y_continuous(
                name = expression(paste(-log[10],"(",italic(P),")")),
                sec.axis = sec_axis(~scale_to_value2(.), name = "Recombination rate (cM/Mb)")) +
            guides(fill = guide_legend(override.aes = list(shape = 22, size = 6))) +
            theme_bw(base_size = 16) +
            theme(
                plot.title = element_text(face = "bold"),
                legend.title.align = 0.5,
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                # legend.margin = margin(1, 1, 1, 1),
                legend.justification = c("right", "top"),
                legend.position = c(0.99, 0.99),
                # legend.spacing = unit(0, "pt"),
                strip.text = element_text(color = "black"),
                strip.text.x = element_blank(),
                # axis.title.y = ggtext::element_markdown(),
                legend.spacing.y = unit(0, "pt")
                ) +
            ggtitle(lead_rsid) +
            xlab(glue::glue("Position on Chromosome {unique(df_locus_r2$CHROM)} (Mb)"))


# Add plot of genes     /work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/test.R

library(AnnotationDbi)
txdb <- loadDb("/work23/home/nsasa/data/CHM13v2.0_hhv6/locusplotr/T2T-CHM13v2_txdb")

txdb
##  TxDb object:
##  # Db type: TxDb
##  # Supporting package: GenomicFeatures
##  # Genome: T2T-CHM13v2.0
##  # Organism: Homo sapiens
##  # Taxonomy ID: 9606
##  # Nb of transcripts: 181708
##  # Db created by: GenomicFeatures package from Bioconductor
##  # Creation time: 2024-02-11 22:40:43 +0900 (Sun, 11 Feb 2024)
##  # GenomicFeatures version at creation time: 1.54.1
##  # RSQLite version at creation time: 2.3.4
##  # DBSCHEMAVERSION: 1.2

# Identify genes in region, create model with exons and gene symbols as labels
gr <- GenomicRanges::GRanges(paste0(lead_chrom), IRanges::IRanges(plot_start, plot_end), strand = "*")
##  GRanges object with 1 range and 0 metadata columns:
##        seqnames            ranges strand
##           <Rle>         <IRanges>  <Rle>
##    [1]       22 50882381-51382381      *
##    -------
##    seqinfo: 1 sequence from an unspecified genome; no seqlengths
gr.txdb <- biovizBase::crunch(txdb, which = gr)
##  Parsing transcripts...
##  Parsing exons...
##  Parsing cds...
##  Parsing utrs...
##  ------exons...
##  ------cdss...
##  ------introns...
##  ------utr...
##  aggregating...
##  Done

gr.txdb
##  GRanges object with 8215 ranges and 4 metadata columns:
##           seqnames            ranges strand |       tx_id     tx_name
##              <Rle>         <IRanges>  <Rle> | <character> <character>
##       [1]       22 50902302-50902351      + |      173204        <NA>
##       [2]       22 50902445-50902687      + |      173204        <NA>
##       [3]       22 50915878-50916064      + |      173204        <NA>
##       [4]       22 50923745-50923882      + |      173204        <NA>
##       [5]       22 50925285-50925350      + |      173204        <NA>
##       ...      ...               ...    ... .         ...         ...
##    [8211]       22 51295792-51295844      * |      175208        <NA>
##    [8212]       22 51296998-51297133      * |      175208        <NA>
##    [8213]       22 51283010-51283372      * |      175209        <NA>
##    [8214]       22 51295792-51295848      * |      175209        <NA>
##    [8215]       22 51296998-51297133      * |      175209        <NA>
##               gene_id     type
##           <character> <factor>
##       [1]      PPP6R2     exon
##       [2]      PPP6R2     exon
##       [3]      PPP6R2     exon
##       [4]      PPP6R2     exon
##       [5]      PPP6R2     exon
##       ...         ...      ...
##    [8211]      RABL2B      utr
##    [8212]      RABL2B      utr
##    [8213]      RABL2B      utr
##    [8214]      RABL2B      utr
##    [8215]      RABL2B      utr
##    -------
##    seqinfo: 1 sequence from T2T-CHM13v2.0 genome

colnames(values(gr.txdb))[4] <- "model"

grl <- split(gr.txdb, gr.txdb$gene_id)
grl <- GenomicRanges::GRangesList(grl, compress = TRUE)

# override angle for plotting gene labels
update_geom_defaults("text", list(angle = 30, hjust = 0))

# plot genes
plot_res <- ggbio::autoplot(grl, aes(type = model), color='brown', fill='brown') +
    theme_light(base_size = 16) +
    scale_y_discrete(expand = expansion(mult = c(0.15, 0.25)))

#   ## fix size (change stepping) by gene_id
#   gene_stepping <- plot_res@ggplot$layers[[3]]$data %>% select(gene_id, stepping) %>% unique()
#   if (length(unique(gene_stepping$stepping))>1) {
#       gene_stepping$stepping_rescaled <- scales::rescale(gene_stepping$stepping, to=c(1,7))
#   } else {
#       gene_stepping$stepping_rescaled <- scales::rescale(gene_stepping$stepping, from=c(0,1), to=c(1,7))
#   }
#   gene_stepping$stepping_rescaled <- scales::rescale(gene_stepping$stepping, to=c(1,7))
#   for (i in 1:8) {
#       for (gene in gene_stepping$gene_id) {
#           plot_res@ggplot$layers[[i]]$data <- plot_res@ggplot$layers[[i]]$data %>%
#                                           mutate(
#                                               stepping=case_when(
#                                                   gene_id == gene ~ gene_stepping[gene_stepping$gene_id==gene,]$stepping_rescaled,
#                                                   TRUE ~ stepping
#                                                   ))
#       }
#   }
#   for (gene in gene_stepping$gene_id) {
#       plot_res@ggplot$layers[[9]]$data <- plot_res@ggplot$layers[[9]]$data %>%
#                                       mutate(
#                                           stepping=case_when(
#                                               .labels == gene ~ gene_stepping[gene_stepping$gene_id==gene,]$stepping_rescaled,
#                                               TRUE ~ stepping
#                                               ))
#   }

#   NumOfLines <- 6
#   max_stepping <- max(plot_res@ggplot$layers[[3]]$data$stepping)
#   if (max_stepping > 1) {
#       for (i in 1:9) {
#           plot_res@ggplot$layers[[i]]$data$stepping <- scales::rescale(plot_res@ggplot$layers[[i]]$data$stepping, from=c(1,max_stepping), to=c(1,1+NumOfLines))
#       }
#   } else {
#       for (i in 1:9) {
#           plot_res@ggplot$layers[[i]]$data$stepping <- scales::rescale(plot_res@ggplot$layers[[i]]$data$stepping, from=c(0,max_stepping), to=c(1,1+NumOfLines))
#       }    
#   }

# return ggplot object with results
gg_geneplot <- plot_res@ggplot +
    theme_light(base_size = 16)
# + coord_cartesian(ylim=c(0.5,7))

geneplot <- gg_geneplot +
    labs(x = glue::glue("Position on Chromosome {lead_chrom} (Mb)")) +
    scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6),
    #limits =c(lead_pos - plot_distance, lead_pos + plot_distance))
    limits =c(min(df_locus_r2$POS), max(df_locus_r2$POS))
    ) +
    theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

regional_assoc_plot_plus <- patchwork::wrap_plots(list(
    regional_assoc_plot +
    labs(x = "") +
    #xlim(lead_pos - plot_distance, lead_pos + plot_distance) +
    xlim(min(df_locus_r2$POS), max(df_locus_r2$POS)) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 0, 5.5)
        ),
    geneplot
    ), nrow = 2, heights = c(2, 1))


#ggsave(regional_assoc_plot_plus, file=paste0(lead_rsid,"_Regional_plot.pdf"), units = "in", height = 8.5, width = 11, device = "pdf")
ggsave(regional_assoc_plot_plus, file=paste0("mafALL", "_", covariate, "_", ADDorDOM, "_each_LD_Regional_plot.pdf"), units = "in", height = 8, width = 12, device = "pdf")








### r2 with eHHV-6B
ehhv6b_ld <- fread(paste0("/work23/home/nsasa/data/CHM13v2.0_hhv6/gwas04_sampleploidy/data04_maf0/eHHV6B_chr22_LDscore.ld"), header=TRUE)
ehhv6b_ld <- ehhv6b_ld[, c("SNP_B", "R2", "DP")]
colnames(ehhv6b_ld) <- c("SNP_B", "R2_eHHV6B", "DP_eHHV6B")
# left_join
df_locus_r2_2 <- left_join(df_locus_r2, ehhv6b_ld, by=c('ID' = 'SNP_B'))
#df_locus_r2_2 <- df_locus_r2_2 %>%
#    mutate(r2_txt = case_when(
#        R2_eHHV6B > 0.2 ~ as.character(paste0("r^2 ~ \"for eHHV-6B = ", sprintf("%0.2f", round(R2_eHHV6B, 2)), "\"")),
#        TRUE ~ NA_character_
#        ))
df_locus_r2_2 <- df_locus_r2_2 %>%
    mutate(r2_shape = case_when(
        R2_eHHV6B > 0.2 ~ "2",
        TRUE ~ "1"
        )) %>%
    mutate(label02 = case_when(
        #R2_eHHV6B > 0.2 ~ as.character(ID),
        P < 1e-4 ~ as.character(ID),
        TRUE ~ NA_character_
        ))

### alpha = sqrt(R2)
regional_assoc_plot_2 <- df_locus_r2_2 %>%
            distinct(ID, .keep_all = TRUE) %>%
            filter(P < plot_pvalue_threshold | R2 > 0.2) %>% # improve pverplotting
            bind_rows(df_locus_r2_2 %>%
                        filter(P >= plot_pvalue_threshold & R2 < 0.2)) %>%
            arrange(desc(color_code)) %>%
            ggplot() +
            geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
            #geom_step(data=recomb_zoom, aes(x=end, y=scale_to_value1(recomb_rate)), color='blue', na.rm=TRUE) +                            # 2nd y axisの最大値をrecomb_rateの最大値に
            geom_step(data=recomb_zoom, aes(x=end, y=scale_to_value1(recomb_rate) * (max(recomb_rate)/100)), color='blue', na.rm=TRUE) +    # 2nd y axisの最大値を100に(locuszoom準拠)
            #geom_point_rast(aes(POS, -log10(P), fill = factor(color_code), size = R2_eHHV6B, alpha = lead, shape = r2_shape)) +
            geom_point_rast(aes(POS, -log10(P), fill = factor(color_code), size = R2_eHHV6B, alpha = sqrt(R2), shape = lead), color='black', stroke=0.5) +
            ggrepel::geom_label_repel(aes(POS, -log10(P), label = label02),
                                        size = 4,
                                        color = "black",
                                        fontface = "bold",
                                        fill = "white",
                                        min.segment.length = 0,
                                        box.padding = 1,
                                        alpha = 1,
                                        nudge_y = 4,
                                        max.overlaps=Inf
                                        ) +
            scale_fill_identity(base::bquote(r^2 ~ " in 1KG"), guide = "legend", labels = c(levels(forcats::fct_drop(df_locus_r2$lead_rsid)), 'No data'), na.translate = FALSE) +
            scale_alpha_identity(base::bquote(r^2 ~ " in 1KG"), guide = "legend", breaks = c(sqrt(0.1), sqrt(0.25), sqrt(0.5), sqrt(1)), labels=c("0.10", "0.25", "0.50", "1.00")) +
            #scale_size_manual(values = c(3, 5), guide = "none") +
            scale_size(name = expression(paste(r^2," for eHHV-6B")),
                breaks = c(0.01, 0.05, 0.1, 0.2, 0.5),
                labels = expression(0.01, 0.05, 0.1, 0.2, 0.5)
                ) +
            scale_shape_manual(values = c(21, 23), guide = "none") +
            #scale_alpha_manual(values = c(0.8, 1), guide = "none") +
            scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            scale_y_continuous(
                name = expression(paste(-log[10],"(",italic(P),")")),
                sec.axis = sec_axis(~scale_to_value2(.), name = "Recombination rate (cM/Mb)")) +
            guides(fill = guide_legend(override.aes = list(shape = 22, size = 6))) +
            theme_bw(base_size = 16) +
            theme(
                plot.title = element_text(face = "bold"),
                legend.title.align = 0.5,
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                # legend.margin = margin(1, 1, 1, 1),
                legend.justification = c("left", "top"),
                legend.position = c(0.01, 0.99),
                legend.direction = "vertical",
                legend.box = "horizontal",
                # legend.spacing = unit(0, "pt"),
                strip.text = element_text(color = "black"),
                strip.text.x = element_blank(),
                # axis.title.y = ggtext::element_markdown(),
                legend.spacing.y = unit(0, "pt"),
                legend.key = element_rect(fill = alpha("white", 0.0))
                ) +
            #ggtitle(lead_rsid) +
            xlab(glue::glue("Position on Chromosome {unique(df_locus_r2$CHROM)} (Mb)"))
            #ggrepel::geom_text_repel(aes(POS, -log10(P), label = r2_txt),
            #                            parse = TRUE,
            #                            size = 4,
            #                            color = "black",
            #                            fontface = "bold",
            #                            min.segment.length = 0,
            #                            box.padding = 1,
            #                            alpha = 1,
            #                            nudge_x = 2
            #                            ) +

regional_assoc_plot_plus_2 <- patchwork::wrap_plots(list(
    regional_assoc_plot_2 +
    labs(x = "") +
    #xlim(lead_pos - plot_distance, lead_pos + plot_distance) +
    xlim(min(df_locus_r2$POS), max(df_locus_r2$POS)) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 0, 5.5)
        ),
    geneplot
    ), nrow = 2, heights = c(2, 1))


#ggsave(regional_assoc_plot_plus, file=paste0(lead_rsid,"_Regional_plot.pdf"), units = "in", height = 8.5, width = 11, device = "pdf")
ggsave(regional_assoc_plot_plus_2, file=paste0("mafALL", "_", covariate, "_", ADDorDOM, "_each_LD_Regional_plot.eHHV6B_LD.pdf"), units = "in", height = 8, width = 12, device = "pdf")

