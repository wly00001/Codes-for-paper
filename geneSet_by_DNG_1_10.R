workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

### distance to neighbor gene (DNG)

library(tidyverse)
library(plyr)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


###input generated tbl
DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))

direction <- c("sd", "su", "od")
plot <- list()
DGE_plotTbl <-list()
DGE_plotTbl_out <-list()
DGE_plotTbl_out2 <-list()

DGE_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
RS_plotTbl <-list()
RS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))
# List.RS <- readRDS("D:/wly/OneDrive - The Wistar Institute/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")

URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
names(URs_clean_genes) <- "gene_symbol"

DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
names(DRas_clean_genes) <- "gene_symbol"

cells <- c("HeLa_1","HeLa_10","HepG2_1","HepG2_10","UP", "U937")
# cells <- c("HeLa_1","HeLa_10","HepG2_1","HepG2_10")

for (d in direction) {
  
  # d="od"
  
  if (d == "sd") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, goi_gene_symbol %in% URs_clean_genes$gene_symbol)
    # DNG <- subset(DNG, nb_gene_symbol %in% DRas_clean_genes$gene_symbol)
  } else if (d == "su") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, nb_gene_symbol %in% URs_clean_genes$gene_symbol)
    # DNG <- subset(DNG, goi_gene_symbol %in% DRas_clean_genes$gene_symbol)
  } else if (d == "od") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, goi_gene_symbol %in% URs_clean_genes$gene_symbol)
  }
  
  KS_test <- list()
  
  for (cell in cells) {
    
    # cell="UP"
    
    DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
    names(DGE) <- c("gene_symbol","L2FC","regu")
    # RS <- List.RS[[cell]][c(1,13)]
    # names(RS) <- c("gene_symbol","deltaRS")
    
    
    df_DGE <- merge(DNG, DGE, by.x="goi_gene_symbol", by.y="gene_symbol")
    
    df_DGE <- df_DGE[order(df_DGE$distance,decreasing=F),]
    bin1=bin2=bin3=bin4=round(nrow(df_DGE)/5)
    bin5=nrow(df_DGE)-(bin1+bin2+bin3+bin4)
    df_DGE$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
    
    DGE_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_DGE, bin == "bin1")),
                                        nrow(subset(df_DGE, bin == "bin2")),
                                        nrow(subset(df_DGE, bin == "bin3")),
                                        nrow(subset(df_DGE, bin == "bin4")),
                                        nrow(subset(df_DGE, bin == "bin5"))),
                             "distance_range"=c(paste(range(subset(df_DGE, bin == "bin1")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin2")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin3")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin4")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin5")$distance),collapse = "~")))
    names(DGE_bin_tbl2) <- paste(names(DGE_bin_tbl2),d,cell,sep = "_")
    DGE_bin_tbl <- cbind(DGE_bin_tbl,DGE_bin_tbl2)
    
    df_DGE3 <- data_summary(df_DGE, varname="L2FC", 
                            groupnames=c("bin"))
    df_DGE3$sampleName <- rep(cell, nrow(df_DGE3))
    DGE_plotTbl[[cell]] <- df_DGE3
    
    df_DGE3_out <- df_DGE3
    df_DGE3_out$sampleName <- rep(paste(d,cell,sep = "_"), nrow(df_DGE3_out))
    DGE_plotTbl_out[[cell]] <- df_DGE3_out
    
    KS_test[[cell]] <- ks.test(subset(df_DGE, bin == "bin1")$L2FC,subset(df_DGE, bin == "bin5")$L2FC)
    
    # df_RS <- merge(DNG, RS, by.x="goi_gene_symbol", by.y="gene_symbol")
    # 
    # df_RS <- df_RS[order(df_RS$distance,decreasing=F),]
    # bin1=bin2=bin3=bin4=round(nrow(df_RS)/5)
    # bin5=nrow(df_RS)-(bin1+bin2+bin3+bin4)
    # df_RS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
    # 
    # RS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_RS, bin == "bin1")),
    #                                    nrow(subset(df_RS, bin == "bin2")),
    #                                    nrow(subset(df_RS, bin == "bin3")),
    #                                    nrow(subset(df_RS, bin == "bin4")),
    #                                    nrow(subset(df_RS, bin == "bin5"))),
    #                         "distance_range"=c(paste(range(subset(df_RS, bin == "bin1")$distance),collapse = "~"),
    #                                            paste(range(subset(df_RS, bin == "bin2")$distance),collapse = "~"),
    #                                            paste(range(subset(df_RS, bin == "bin3")$distance),collapse = "~"),
    #                                            paste(range(subset(df_RS, bin == "bin4")$distance),collapse = "~"),
    #                                            paste(range(subset(df_RS, bin == "bin5")$distance),collapse = "~")))
    # names(RS_bin_tbl2) <- paste(names(RS_bin_tbl2),d,cell,sep = "_")
    # RS_bin_tbl <- cbind(RS_bin_tbl,RS_bin_tbl2)
    # 
    # df_RS3 <- data_summary(df_RS, varname="deltaRS", 
    #                        groupnames=c("bin"))
    # df_RS3$sampleName <- rep(cell, nrow(df_RS3))
    # RS_plotTbl[[cell]] <- df_RS3
    
    
  }
  
  
  DGE_plotTbl_JTE <- do.call(rbind,DGE_plotTbl)
  DGE_plotTbl_out2[[d]] <- do.call(rbind,DGE_plotTbl_out)
  
  library(RColorBrewer)
  plot[[paste0(d,"_DGE")]] <- ggplot(DGE_plotTbl_JTE, aes(x=bin, y=L2FC, group=sampleName, color=sampleName)) + 
    geom_line(size=1) +
    geom_point(size=3) +
    # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
    geom_errorbar(aes(ymin=L2FC-se, ymax=L2FC+se), width=0.5, size=1,
                  position=position_dodge(0.05)) +
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    labs(title = paste0(d,
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(1)],"=",round(KS_test[[1]]$statistic,2),", pval=",formatC(KS_test[[1]]$p.value, format = "e", digits = 2),
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(2)],"=",round(KS_test[[2]]$statistic,2),", pval=",formatC(KS_test[[2]]$p.value, format = "e", digits = 2),
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(3)],"=",round(KS_test[[3]]$statistic,2),", pval=",formatC(KS_test[[3]]$p.value, format = "e", digits = 2),
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(4)],"=",round(KS_test[[4]]$statistic,2),", pval=",formatC(KS_test[[4]]$p.value, format = "e", digits = 2),
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(5)],"=",round(KS_test[[5]]$statistic,2),", pval=",formatC(KS_test[[5]]$p.value, format = "e", digits = 2),
                        "\nKS_bin1.v.bin5_",names(KS_test)[c(6)],"=",round(KS_test[[6]]$statistic,2),", pval=",formatC(KS_test[[6]]$p.value, format = "e", digits = 2)),
         x = "distance size bin", y = "Mean of L2FC") +
    scale_color_manual(values=c(#"gray80","gray50","gray0",
      brewer.pal(n = 12,name = 'Paired')[c(1,2,3,4,5,7)]))
  
  
  # RS_plotTbl_JTE <- do.call(rbind,RS_plotTbl)
  # 
  # library(RColorBrewer)
  # plot[[paste0(d,"_RS")]] <- ggplot(RS_plotTbl_JTE, aes(x=bin, y=deltaRS, group=sampleName, color=sampleName)) + 
  #   geom_line(size=1) +
  #   geom_point(size=3) +
  #   # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
  #   geom_errorbar(aes(ymin=deltaRS-se, ymax=deltaRS+se), width=0.5, size=1,
  #                 position=position_dodge(0.05)) +
  #   theme_classic()+
  #   theme(axis.text.y=element_text(colour = "black"),
  #         axis.title.y = element_text(colour = "black"),
  #         axis.text.x=element_text(colour = "black"),
  #         axis.title.x = element_text(colour = "black"),
  #         axis.line = element_line(colour = "black"),
  #         axis.ticks = element_line(size = 1),
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  #   labs(title = paste0(d),
  #        x = "distance size bin", y = "Mean of deltaRS") +
  #   scale_color_manual(values=c(#"gray80","gray50","gray0",
  #     brewer.pal(n = 12,name = 'Paired')[c(1,3,5,7)]))
  
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 3)

# write.csv(as.data.frame(t(DGE_bin_tbl)), "./_bin_table/DGE_bin_tbl.csv", col.names = F, row.names = T)
# write.csv(as.data.frame(t(RS_bin_tbl)), "./_bin_table/RS_bin_tbl.csv", col.names = F, row.names = T)

DGE_plotTbl_out3 <- do.call(rbind,DGE_plotTbl_out2)
# write.csv(DGE_plotTbl_out3, "./_bin_table/DGE_distance_bin_L2FC_tbl.csv", col.names = F, row.names = T)








############################################################################################################
############################################################################################################
############################################################################################################

### Gene expression grid analysis ~ L2FC, goi vs NG (opposite downstream neighbor gene)

library(tidyverse)
library(plyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))


# DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/DNG/DNGlist.RData"))
DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))

# direction <- c("sd", "su", "od")

direction <- c("od")
dir.create(file.path("./", "_grid"))
hm_plot <- list()
DGE_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
library(openxlsx)
wb = createWorkbook()

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))
# List.RS <- readRDS("D:/wly/OneDrive - The Wistar Institute/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")

URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
names(URs_clean_genes) <- "gene_symbol"

DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
names(DRas_clean_genes) <- "gene_symbol"

cells <- c("HeLa_1","HeLa_10","HepG2_1","HepG2_10","U937","UP")

for (d in direction) {
  
  # d="su"
  
  if (d == "sd") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, goi_gene_symbol %in% URs_clean_genes$gene_symbol)
    # DNG <- subset(DNG, nb_gene_symbol %in% DRas_clean_genes$gene_symbol)
  } else if (d == "su") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, nb_gene_symbol %in% URs_clean_genes$gene_symbol)
    # DNG <- subset(DNG, goi_gene_symbol %in% DRas_clean_genes$gene_symbol)
  } else if (d == "od") {
    DNG <- DNGlist[[d]]
    # DNG <- subset(DNG, goi_gene_symbol %in% URs_clean_genes$gene_symbol)
  }
  
  ALL_grid_bin_tbl <- list()
  
  for (cell in cells) {
    
    # cell="HepG2_10"
    grid_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
    
    DGE_goi <- List.DGE_QSF[[cell]][c(1,4,6,10)]
    names(DGE_goi) <- c("gene_symbol_goi","rpm_DMSO_goi","L2FC_goi","regu_goi")
    DGE_NG <- List.DGE_QSF[[cell]][c(1,4,6,10)]
    names(DGE_NG) <- c("gene_symbol_NG","rpm_DMSO_NG","L2FC_NG","regu_NG")
    
    
    df_DGE <- merge(DNG, DGE_goi, by.x="goi_gene_symbol", by.y="gene_symbol_goi")
    df_DGE <- merge(df_DGE, DGE_NG, by.x="nb_gene_symbol", by.y="gene_symbol_NG")
    
    
    df_DGE <- df_DGE[order(df_DGE$distance,decreasing=F),]
    bin1=bin2=bin3=bin4=round(nrow(df_DGE)/5)
    bin5=nrow(df_DGE)-(bin1+bin2+bin3+bin4)
    df_DGE$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
    
    DGE_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_DGE, bin == "bin1")),
                                        nrow(subset(df_DGE, bin == "bin2")),
                                        nrow(subset(df_DGE, bin == "bin3")),
                                        nrow(subset(df_DGE, bin == "bin4")),
                                        nrow(subset(df_DGE, bin == "bin5"))),
                             "distance_range"=c(paste(range(subset(df_DGE, bin == "bin1")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin2")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin3")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin4")$distance),collapse = "~"),
                                                paste(range(subset(df_DGE, bin == "bin5")$distance),collapse = "~")))
    names(DGE_bin_tbl2) <- paste(names(DGE_bin_tbl2),d,cell,sep = "_")
    DGE_bin_tbl <- cbind(DGE_bin_tbl,DGE_bin_tbl2)
    
    
    
    b1="bin1"
    b2="bin2"
    
    
    plotDF <- subset(df_DGE, bin==b1 | bin==b2)
    
    df2D <- plotDF
    bin1=bin2=bin3=bin4=round(nrow(df2D)/5)
    bin5=nrow(df2D)-(bin1+bin2+bin3+bin4)
    
    df2D <- df2D[order(df2D$rpm_DMSO_goi,decreasing=F),]
    df2D$BIN_rpm_DMSO_goi <- c(rep("Bin1",bin1),rep("Bin2",bin2),rep("Bin3",bin3),rep("Bin4",bin4),rep("Bin5",bin5))
    
    df2D <- df2D[order(df2D$rpm_DMSO_NG,decreasing=F),]
    df2D$BIN_rpm_DMSO_NG <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
    
    
    df2D$bin <- paste0(df2D$BIN_rpm_DMSO_goi,"_",df2D$BIN_rpm_DMSO_NG)
    
    range_bin1bin2 <- df2D %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_goi=length(rpm_DMSO_goi),
                       rpm_min_goi=min(rpm_DMSO_goi),
                       rpm_max_goi=max(rpm_DMSO_goi),
                       gene_N_NG=length(rpm_DMSO_NG),
                       rpm_min_NG=min(rpm_DMSO_NG),
                       rpm_max_NG=max(rpm_DMSO_NG))
    
    range_bin1bin2$goi_rpm_range <- paste(round(range_bin1bin2$rpm_min_goi,2),round(range_bin1bin2$rpm_max_goi,2),sep = "-")
    range_bin1bin2$NG_rpm_range <- paste(round(range_bin1bin2$rpm_min_NG,2),round(range_bin1bin2$rpm_max_NG,2),sep = "-")
    range_bin1bin2$gene_N <- range_bin1bin2$gene_N_goi
    range_bin1bin2$gene_N.goi_rpm_range.NG_rpm_range <- paste(range_bin1bin2$gene_N,range_bin1bin2$goi_rpm_range,range_bin1bin2$NG_rpm_range,sep = ";")
    range_bin1bin2 <- range_bin1bin2[c(1,11)] %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    range_bin1bin2 <- dcast(range_bin1bin2, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "gene_N.goi_rpm_range.NG_rpm_range")
    range_bin1bin2 <- range_bin1bin2[order(nrow(range_bin1bin2):1),]
    
    
    df2D0 <- df2D
    
    MedianAll <- median(df2D0$L2FC_goi)
    df2D2 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(Median=median(L2FC_goi))
    df2D2$value <- df2D2$Median #- MedianAll
    df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    df2D3 <- dcast(df2D2, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_rpm_DMSO_goi)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot[[paste0(d,".",cell,".",b1,b2,".goi")]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.L2FC", column_title = paste0(d,".",cell,".",b1,b2,".goi"), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    # write.csv(matrix2D, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig4i.csv"),row.names = T)
    
    
    
    MedianAll <- median(df2D0$L2FC_NG)
    df2D2 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(Median=median(L2FC_NG))
    df2D2$value <- df2D2$Median #- MedianAll
    df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    df2D3 <- dcast(df2D2, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_rpm_DMSO_goi)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot[[paste0(d,".",cell,".",b1,b2,".NG")]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.L2FC", column_title = paste0(d,".",cell,".",b1,b2,".NG"), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
    
    
    
    b1="bin4"
    b2="bin5"
    
    
    plotDF <- subset(df_DGE, bin==b1 | bin==b2)
    
    df2D <- plotDF
    bin1=bin2=bin3=bin4=round(nrow(df2D)/5)
    bin5=nrow(df2D)-(bin1+bin2+bin3+bin4)
    
    df2D <- df2D[order(df2D$rpm_DMSO_goi,decreasing=F),]
    df2D$BIN_rpm_DMSO_goi <- c(rep("Bin1",bin1),rep("Bin2",bin2),rep("Bin3",bin3),rep("Bin4",bin4),rep("Bin5",bin5))
    
    df2D <- df2D[order(df2D$rpm_DMSO_NG,decreasing=F),]
    df2D$BIN_rpm_DMSO_NG <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
    
    
    df2D$bin <- paste0(df2D$BIN_rpm_DMSO_goi,"_",df2D$BIN_rpm_DMSO_NG)
    
    range_bin4bin5 <- df2D %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_goi=length(rpm_DMSO_goi),
                       rpm_min_goi=min(rpm_DMSO_goi),
                       rpm_max_goi=max(rpm_DMSO_goi),
                       gene_N_NG=length(rpm_DMSO_NG),
                       rpm_min_NG=min(rpm_DMSO_NG),
                       rpm_max_NG=max(rpm_DMSO_NG))
    
    range_bin4bin5$goi_rpm_range <- paste(round(range_bin4bin5$rpm_min_goi,2),round(range_bin4bin5$rpm_max_goi,2),sep = "-")
    range_bin4bin5$NG_rpm_range <- paste(round(range_bin4bin5$rpm_min_NG,2),round(range_bin4bin5$rpm_max_NG,2),sep = "-")
    range_bin4bin5$gene_N <- range_bin4bin5$gene_N_goi
    range_bin4bin5$gene_N.goi_rpm_range.NG_rpm_range <- paste(range_bin4bin5$gene_N,range_bin4bin5$goi_rpm_range,range_bin4bin5$NG_rpm_range,sep = ";")
    range_bin4bin5 <- range_bin4bin5[c(1,11)] %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    range_bin4bin5 <- dcast(range_bin4bin5, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "gene_N.goi_rpm_range.NG_rpm_range")
    range_bin4bin5 <- range_bin4bin5[order(nrow(range_bin4bin5):1),]
    
    ALL_grid_bin_tbl[[cell]] <- cbind(range_bin1bin2,range_bin4bin5)
    
    
    df2D0 <- df2D
    
    MedianAll <- median(df2D0$L2FC_goi)
    df2D2 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(Median=median(L2FC_goi))
    df2D2$value <- df2D2$Median #- MedianAll
    df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    df2D3 <- dcast(df2D2, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_rpm_DMSO_goi)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot[[paste0(d,".",cell,".",b1,b2,".goi")]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.L2FC", column_title = paste0(d,".",cell,".",b1,b2,".goi"), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
    
    MedianAll <- median(df2D0$L2FC_NG)
    df2D2 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(Median=median(L2FC_NG))
    df2D2$value <- df2D2$Median #- MedianAll
    df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_rpm_DMSO_goi', 'BIN_rpm_DMSO_NG'), remove = F)
    df2D3 <- dcast(df2D2, BIN_rpm_DMSO_goi~BIN_rpm_DMSO_NG, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_rpm_DMSO_goi)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot[[paste0(d,".",cell,".",b1,b2,".NG")]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.L2FC", column_title = paste0(d,".",cell,".",b1,b2,".NG"), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
    
    
  }
  
  DF_ALL_grid_bin_tbl <- do.call(rbind,ALL_grid_bin_tbl)
  
  # write.csv(as.data.frame(DF_ALL_grid_bin_tbl), paste0("./_bin_table/GEgrid_GE_bin_tbl_",d,".csv"), col.names = F, row.names = T)
  
  name <- d
  
  # Write to Excel
  addWorksheet(wb, name)
  writeData(wb, sheet = d, as.data.frame(DF_ALL_grid_bin_tbl), rowNames = T, colNames = T)
  
  
}

# write.csv(as.data.frame(t(DGE_bin_tbl)), "./_bin_table/GEgrid_geneDistance_bin_tbl.csv", col.names = F, row.names = T)
# saveWorkbook(wb, "./_bin_table/GEgrid_GE_bin_tbl.xlsx", overwrite = T)


library(cowplot)
plot_grid(plotlist=hm_plot[c(1:16)],ncol = 4)
# plot_grid(plotlist=hm_plot[c(17:32)],ncol = 4)
# plot_grid(plotlist=hm_plot[c(33:48)],ncol = 4)
# plot_grid(plotlist=hm_plot[c(49:64)],ncol = 4)
plot_grid(plotlist=hm_plot[c(9)],ncol = 1)
plot_grid(plotlist=hm_plot[c(13)],ncol = 1)
plot_grid(plotlist=hm_plot[c(17)],ncol = 1)
plot_grid(plotlist=hm_plot[c(21)],ncol = 1)









