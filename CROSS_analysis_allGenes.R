workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

library(tidyverse)
library(RColorBrewer)

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CROSS1 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS1_1000bp.RData"))
List.CROSS2 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS2.RData"))

List.CROSSa <- List.CROSS1
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS1","CROSS",names(x)); x })
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS2","CROSS",names(x)); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CROSSa <- lapply(List.CROSSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })

ggplot() +
  stat_ecdf(data=List.CROSSa[[c(1)]], aes(x=CROSS.JTE, color = "a"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSa[[c(2)]], aes(x=CROSS.JTE, color = "b"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSa[[c(3)]], aes(x=CROSS.JTE, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSa[[c(4)]], aes(x=CROSS.JTE, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
    brewer.pal(n = 12,name = 'Paired')),
    labels = c(paste0(names(List.CROSSa)[c(1)],", ",nrow(List.CROSSa[[c(1)]]),", ",round(median(List.CROSSa[[c(1)]][[c(4)]]),2)),
               paste0(names(List.CROSSa)[c(2)],", ",nrow(List.CROSSa[[c(2)]]),", ",round(median(List.CROSSa[[c(2)]][[c(4)]]),2)),
               paste0(names(List.CROSSa)[c(3)],", ",nrow(List.CROSSa[[c(3)]]),", ",round(median(List.CROSSa[[c(3)]][[c(4)]]),2)),
               paste0(names(List.CROSSa)[c(4)],", ",nrow(List.CROSSa[[c(4)]]),", ",round(median(List.CROSSa[[c(4)]][[c(4)]]),2)))) +
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        plot.title = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(xlim = c(-8, 2)) +
  labs(title = paste("CROSS, all UR-isolated genes",sep = ""),
       x = "CROSS", y = "Cumulative fraction", color = "")





List.DGE_QSF2 <- lapply(List.DGE_QSF, function(x) { x <- x[c(1,4,5)]; x })
DGE_QSF <- Reduce(function(x, y) merge(x, y), List.DGE_QSF2)
DGE_QSF2 <- subset(DGE_QSF, DGE_QSF[[2]]>=10 & DGE_QSF[[3]]>=10 & DGE_QSF[[4]]>=10 & DGE_QSF[[5]]>=10 & DGE_QSF[[6]]>=10 & DGE_QSF[[7]]>=10 & DGE_QSF[[8]]>=10 & DGE_QSF[[9]]>=10)


List.CROSSb <- lapply(List.CROSSa, function(x) { x <- subset(x,gene_symbol %in% DGE_QSF2$gene_symbol); x })

ggplot() +
  stat_ecdf(data=List.CROSSb[[c(1)]], aes(x=CROSS.JTE, color = "a"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSb[[c(2)]], aes(x=CROSS.JTE, color = "b"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSb[[c(3)]], aes(x=CROSS.JTE, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=List.CROSSb[[c(4)]], aes(x=CROSS.JTE, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
    brewer.pal(n = 12,name = 'Paired')),
    labels = c(paste0(names(List.CROSSb)[c(1)],", ",nrow(List.CROSSb[[c(1)]]),", ",round(median(List.CROSSb[[c(1)]][[c(4)]]),2)),
               paste0(names(List.CROSSb)[c(2)],", ",nrow(List.CROSSb[[c(2)]]),", ",round(median(List.CROSSb[[c(2)]][[c(4)]]),2)),
               paste0(names(List.CROSSb)[c(3)],", ",nrow(List.CROSSb[[c(3)]]),", ",round(median(List.CROSSb[[c(3)]][[c(4)]]),2)),
               paste0(names(List.CROSSb)[c(4)],", ",nrow(List.CROSSb[[c(4)]]),", ",round(median(List.CROSSb[[c(4)]][[c(4)]]),2)))) +
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        plot.title = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(xlim = c(-8, 2)) +
  labs(title = paste("CROSS, commonly expressed UR-isolated genes (RPM >=10 in all QSF samples)",sep = ""),
       x = "CROSS", y = "Cumulative fraction", color = "")




#### QSF_DGE vs readthru


library(tidyverse)
library(RColorBrewer)
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


List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CROSS1 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS1_1000bp.RData"))
List.CROSS2 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS2.RData"))

List.CROSSa <- List.CROSS1
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS1","CROSS",names(x)); x })
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS2","CROSS",names(x)); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CROSSa <- lapply(List.CROSSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
scatterplot <- list()

CROSS_plotTbl <-list()
CROSS_plotTbl_out <-list()
CROSS_plotTbl_out2 <-list()
CROSS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
KS_test <- list()


for (cell in cells) {
  
  # cell="HepG2"
  
  ldf <- List.DGE_QSF[[cell]]
  QSF_UP <- subset(ldf,ldf[[10]]=="UP")
  QSF_NC <- subset(ldf,ldf[[10]]=="NO")
  QSF_DN <- subset(ldf,ldf[[10]]=="DN")
  
  
  CROSS <- List.CROSSa[[cell]]
  
  CDFplot[[cell]] <- ggplot() +
    stat_ecdf(data=subset(CROSS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=CROSS.JTE, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(CROSS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=CROSS.JTE, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(CROSS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=CROSS.JTE, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(cell, ", QSF.DGE vs CROSS",
                        "\nNCvUP, ", formatC(ks.test(subset(CROSS, gene_symbol %in% QSF_NC$gene_symbol)$CROSS.JTE, subset(CROSS, gene_symbol %in% QSF_UP$gene_symbol)$CROSS.JTE)$p.value, format = "e", digits = 2),
                        "\nNCvDN, ", formatC(ks.test(subset(CROSS, gene_symbol %in% QSF_NC$gene_symbol)$CROSS.JTE, subset(CROSS, gene_symbol %in% QSF_DN$gene_symbol)$CROSS.JTE)$p.value, format = "e", digits = 2)),
         x = "CROSS", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(CROSS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(CROSS, gene_symbol %in% QSF_NC$gene_symbol)$CROSS.JTE),2)), 
                                  paste0("UP genes, ",nrow(subset(CROSS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(CROSS, gene_symbol %in% QSF_UP$gene_symbol)$CROSS.JTE),2)),
                                  paste0("DN genes, ",nrow(subset(CROSS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(CROSS, gene_symbol %in% QSF_DN$gene_symbol)$CROSS.JTE),2))
                                  
    ), 
    values=c("grey","red","blue")) + 
    coord_cartesian(xlim = c(-8, 2))
  
  
  DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
  names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
  CROSS <- List.CROSSa[[cell]]
  DF <- merge(DGE,CROSS)
  
  cor_all <- cor.test(DF$log2Ratio,DF$CROSS.JTE)
  scatterplot[[cell]] <- ggplot(data=DF, aes(x=log2Ratio, y=CROSS.JTE, label=gene_symbol)) +
    geom_point(size = 1.5, alpha=0.3) +
    # geom_bin2d(bins = 300) +
    scale_fill_gradient(low="darkblue",high="red") +
    theme_bw() + 
    # coord_fixed(ratio = 1)+
    # scale_x_continuous(limits = c(-3, 3))+
    # scale_y_continuous(limits = c(-3, 3))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(title = paste0(cell, ", cor.all=",round(cor_all$estimate,2),", pval=",formatC(cor_all$p.value, format = "e", digits = 2)),
         x = "log2Ratio", y = "CROSS.JTE")
  
  
  
  
  
  
  
  ######### 5 bins by GOI expression level in DMSO 
  
  DGE <- List.DGE_QSF[[cell]][c(1,4,6,10)]
  names(DGE) <- c("gene_symbol","RPM_DMSO","log2Ratio","regu.adj")
  CROSS <- List.CROSSa[[cell]]
  df_CROSS <- merge(DGE,CROSS)
  
  df_CROSS <- df_CROSS[order(df_CROSS$RPM_DMSO,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CROSS)/5)
  bin5=nrow(df_CROSS)-(bin1+bin2+bin3+bin4)
  df_CROSS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CROSS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CROSS, bin == "bin1")),
                                        nrow(subset(df_CROSS, bin == "bin2")),
                                        nrow(subset(df_CROSS, bin == "bin3")),
                                        nrow(subset(df_CROSS, bin == "bin4")),
                                        nrow(subset(df_CROSS, bin == "bin5"))),
                             "RPM_DMSO_range"=c(paste(range(subset(df_CROSS, bin == "bin1")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin2")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin3")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin4")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin5")$RPM_DMSO),collapse = "~")))
  names(CROSS_bin_tbl2) <- paste(names(CROSS_bin_tbl2),cell,sep = "_")
  CROSS_bin_tbl <- cbind(CROSS_bin_tbl,CROSS_bin_tbl2)
  
  for (b in c("bin1","bin2","bin3","bin4","bin5")){
    
    # b="bin1"
    
    ldf <- subset(df_CROSS,bin==b)
    
    QSF_UP <- subset(ldf,regu.adj=="UP")
    QSF_NC <- subset(ldf,regu.adj=="NO")
    QSF_DN <- subset(ldf,regu.adj=="DN")
    
    
    CDFplot[[paste0(cell,"_RPM_DMSO_",b)]] <- ggplot() +
      stat_ecdf(data=QSF_NC, aes(x=CROSS.JTE, color = "a"),geom = "step", size=1) +
      stat_ecdf(data=QSF_UP, aes(x=CROSS.JTE, color = "b"),geom = "step", size=1) +
      stat_ecdf(data=QSF_DN, aes(x=CROSS.JTE, color = "c"),geom = "step", size=1) +
      theme_bw() +
      theme(plot.title = element_text(size = 10)) +
      labs(title = paste0(cell,"_RPM_DMSO_",b, ", QSF.DGE vs CROSS",
                          "\nNCvUP, ", formatC(ks.test(QSF_NC$CROSS.JTE, QSF_UP$CROSS.JTE)$p.value, format = "e", digits = 2),
                          "\nNCvDN, ", formatC(ks.test(QSF_NC$CROSS.JTE, QSF_DN$CROSS.JTE)$p.value, format = "e", digits = 2)),
           x = "CROSS", y = "Cumulative fraction", color = "") +
      scale_color_manual(labels = c(paste0("NC genes, ",nrow(QSF_NC),", ",round(median(QSF_NC$CROSS.JTE),2)), 
                                    paste0("UP genes, ",nrow(QSF_UP),", ",round(median(QSF_UP$CROSS.JTE),2)),
                                    paste0("DN genes, ",nrow(QSF_DN),", ",round(median(QSF_DN$CROSS.JTE),2))
                                    
      ), 
      values=c("grey","red","blue")) + 
      coord_cartesian(xlim = c(-8, 2))
    
    
  }
}

library(cowplot)
plot_grid(plotlist=CDFplot,ncol = 6)
plot_grid(plotlist=CDFplot[c(19)])
plot_grid(plotlist=CDFplot[c(20:24)],ncol = 3)

# plot_grid(plotlist=scatterplot,ncol = 4)

# write.csv(merge(List.CROSSa[["HepG2"]][c(1,4)],List.DGE_QSF[["HepG2"]][c(1,10)]), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig4d.csv"),row.names = T)






#### CROSS vs distance

library(tidyverse)
library(RColorBrewer)
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

DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))
DNG_od <- DNGlist[["od"]]

plot <- list()
CROSS_plotTbl <-list()
CROSS_plotTbl_out <-list()
CROSS_plotTbl_out2 <-list()

CROSS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

DNG <- DNG_od

List.CROSS1 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS1_1000bp.RData"))
List.CROSS2 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS2.RData"))

List.CROSSa <- List.CROSS1
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS1","CROSS",names(x)); x })
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS2","CROSS",names(x)); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CROSSa <- lapply(List.CROSSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })



cells <- c("UP", "U937", "HeLa", "HepG2")

KS_test <- list()

for (cell in cells) {
  
  # cell="UP"
  
  CROSS <- List.CROSSa[[cell]][c(1,4)]
  names(CROSS) <- c("gene_symbol","CROSS")
  
  
  df_CROSS <- merge(DNG, CROSS, by.x="goi_gene_symbol", by.y="gene_symbol")
  
  df_CROSS <- df_CROSS[order(df_CROSS$distance,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CROSS)/5)
  bin5=nrow(df_CROSS)-(bin1+bin2+bin3+bin4)
  df_CROSS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CROSS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CROSS, bin == "bin1")),
                                        nrow(subset(df_CROSS, bin == "bin2")),
                                        nrow(subset(df_CROSS, bin == "bin3")),
                                        nrow(subset(df_CROSS, bin == "bin4")),
                                        nrow(subset(df_CROSS, bin == "bin5"))),
                             "distance_range"=c(paste(range(subset(df_CROSS, bin == "bin1")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin2")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin3")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin4")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin5")$distance),collapse = "~")))
  names(CROSS_bin_tbl2) <- paste(names(CROSS_bin_tbl2),cell,sep = "_")
  CROSS_bin_tbl <- cbind(CROSS_bin_tbl,CROSS_bin_tbl2)
  
  df_CROSS3 <- data_summary(df_CROSS, varname="CROSS", 
                            groupnames=c("bin"))
  df_CROSS3$sampleName <- rep(cell, nrow(df_CROSS3))
  CROSS_plotTbl[[cell]] <- df_CROSS3
  
  df_CROSS3_out <- df_CROSS3
  df_CROSS3_out$sampleName <- rep(paste(cell,sep = "_"), nrow(df_CROSS3_out))
  CROSS_plotTbl_out[[cell]] <- df_CROSS3_out
  
  KS_test[[cell]] <- ks.test(subset(df_CROSS, bin == "bin1")$CROSS,subset(df_CROSS, bin == "bin5")$CROSS)
  
  
  
}


CROSS_plotTbl_JTE <- do.call(rbind,CROSS_plotTbl)
CROSS_plotTbl_out2 <- do.call(rbind,CROSS_plotTbl_out)
# write.csv(CROSS_plotTbl_out2, "./_bin_table/CROSS_plotTbl.csv", col.names = F, row.names = T)

library(RColorBrewer)
ggplot(CROSS_plotTbl_JTE, aes(x=bin, y=CROSS, group=sampleName, color=sampleName)) + 
  geom_line(size=1) +
  geom_point(size=3) +
  # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
  geom_errorbar(aes(ymin=CROSS-se, ymax=CROSS+se), width=0.5, size=1,
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
  labs(title = paste0("KS_bin1.v.bin5_",names(KS_test)[c(1)],"=",round(KS_test[[1]]$statistic,2),", pval=",formatC(KS_test[[1]]$p.value, format = "e", digits = 2),
                      "\nKS_bin1.v.bin5_",names(KS_test)[c(2)],"=",round(KS_test[[2]]$statistic,2),", pval=",formatC(KS_test[[2]]$p.value, format = "e", digits = 2),
                      "\nKS_bin1.v.bin5_",names(KS_test)[c(3)],"=",round(KS_test[[3]]$statistic,2),", pval=",formatC(KS_test[[3]]$p.value, format = "e", digits = 2),
                      "\nKS_bin1.v.bin5_",names(KS_test)[c(4)],"=",round(KS_test[[4]]$statistic,2),", pval=",formatC(KS_test[[4]]$p.value, format = "e", digits = 2)),
       x = "distance size bin", y = "Mean of CROSS") +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
    brewer.pal(n = 12,name = 'Paired')[c(1,2,3,4)]))











#### DGE vs CROSS by distance

library(tidyverse)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

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

DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))
DNG_od <- DNGlist[["od"]]

plot <- list()
CROSS_plotTbl <-list()
CROSS_plotTbl_out <-list()
CROSS_plotTbl_out2 <-list()

CROSS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

DNG <- DNG_od

List.CROSS1 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS1_1000bp.RData"))
List.CROSS2 <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CROSS_table/List.CROSS2.RData"))

List.CROSSa <- List.CROSS1
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS1","CROSS",names(x)); x })
List.CROSSa <- lapply(List.CROSSa, function(x) { names(x) <- gsub("CROSS2","CROSS",names(x)); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CROSSa <- lapply(List.CROSSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })


List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")


cells <- c("UP", "U937", "HeLa", "HepG2")

KS_test <- list()
scatterplot <- list()
hm_plot <- list()
hm_plot_all <- list()

for (cell in cells) {
  
  # cell="HepG2"
  
  CROSS <- List.CROSSa[[cell]][c(1,4)]
  names(CROSS) <- c("gene_symbol","CROSS")
  
  
  df_CROSS <- merge(DNG, CROSS, by.x="goi_gene_symbol", by.y="gene_symbol")
  
  df_CROSS <- df_CROSS[order(df_CROSS$distance,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CROSS)/5)
  bin5=nrow(df_CROSS)-(bin1+bin2+bin3+bin4)
  df_CROSS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CROSS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CROSS, bin == "bin1")),
                                        nrow(subset(df_CROSS, bin == "bin2")),
                                        nrow(subset(df_CROSS, bin == "bin3")),
                                        nrow(subset(df_CROSS, bin == "bin4")),
                                        nrow(subset(df_CROSS, bin == "bin5"))),
                             "distance_range"=c(paste(range(subset(df_CROSS, bin == "bin1")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin2")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin3")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin4")$distance),collapse = "~"),
                                                paste(range(subset(df_CROSS, bin == "bin5")$distance),collapse = "~")))
  names(CROSS_bin_tbl2) <- paste(names(CROSS_bin_tbl2),cell,sep = "_")
  CROSS_bin_tbl <- cbind(CROSS_bin_tbl,CROSS_bin_tbl2)
  
  df_CROSS3 <- data_summary(df_CROSS, varname="CROSS", 
                            groupnames=c("bin"))
  df_CROSS3$sampleName <- rep(cell, nrow(df_CROSS3))
  CROSS_plotTbl[[cell]] <- df_CROSS3
  
  df_CROSS3_out <- df_CROSS3
  df_CROSS3_out$sampleName <- rep(paste(cell,sep = "_"), nrow(df_CROSS3_out))
  CROSS_plotTbl_out[[cell]] <- df_CROSS3_out
  
  KS_test[[cell]] <- ks.test(subset(df_CROSS, bin == "bin1")$CROSS,subset(df_CROSS, bin == "bin5")$CROSS)
  
  
  for (b in c("bin1","bin2","bin3","bin4","bin5")){
    
    # b="bin1"
    
    DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
    names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
    CROSS_bin <- subset(df_CROSS,bin==b)[c(1,11,12)]
    names(CROSS_bin) <- c("gene_symbol","distance","CROSS.JTE")
    DF <- merge(DGE,CROSS_bin)
    
    cor_all <- cor.test(DF$log2Ratio,DF$CROSS.JTE)
    scatterplot[[paste0(cell,"_",b)]] <- ggplot(data=DF, aes(x=log2Ratio, y=CROSS.JTE, label=gene_symbol)) +
      geom_point(size = 1.5, alpha=0.3) +
      # geom_bin2d(bins = 300) +
      scale_fill_gradient(low="darkblue",high="red") +
      theme_bw() + 
      theme(plot.title = element_text(size=12)) +
      # coord_fixed(ratio = 1)+
      # scale_x_continuous(limits = c(-3, 3))+
      # scale_y_continuous(limits = c(-3, 3))+
      geom_vline(xintercept=0, col="black", linetype="solid") +
      geom_hline(yintercept=0, col="black", linetype="solid") +
      labs(title = paste0(cell,", ",b, ", cor.all=",round(cor_all$estimate,2),", pval=",formatC(cor_all$p.value, format = "e", digits = 2)),
           x = "log2Ratio", y = "CROSS.JTE")
    
    
    
    
    
    
    
    
    
    df2D <- DF
    binA1=binA2=binA3=binA4=round(nrow(df2D)/5)
    binA5=nrow(df2D)-(binA1+binA2+binA3+binA4)
    
    df2D <- df2D[order(df2D$distance,decreasing=F),]
    df2D$BIN_distance <- c(rep("BinD1",binA1),rep("BinD2",binA2),rep("BinD3",binA3),rep("BinD4",binA4),rep("BinD5",binA5))
    
    df2D <- df2D[order(df2D$CROSS.JTE,decreasing=F),]
    df2D$BIN_CROSS.JTE <- c(rep("binC1",binA1),rep("binC2",binA2),rep("binC3",binA3),rep("binC4",binA4),rep("binC5",binA5))
    
    
    df2D$bin <- paste0(df2D$BIN_distance,"_",df2D$BIN_CROSS.JTE)
    
    range_bin <- df2D %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_distance=length(distance),
                       min_distance=min(distance),
                       max_distance=max(distance),
                       gene_N_CROSS.JTE=length(CROSS.JTE),
                       min_CROSS.JTE=min(CROSS.JTE),
                       max_CROSS.JTE=max(CROSS.JTE))
    
    range_bin$distance_range <- paste(round(range_bin$min_distance,2),round(range_bin$max_distance,2),sep = " ~ ")
    range_bin$CROSS.JTE_range <- paste(round(range_bin$min_CROSS.JTE,2),round(range_bin$max_CROSS.JTE,2),sep = " ~ ")
    range_bin$gene_N <- range_bin$gene_N_distance
    range_bin$gene_N.distance_range.CROSS.JTE_range <- paste(range_bin$gene_N,range_bin$distance_range,range_bin$CROSS.JTE_range,sep = ";")
    range_bin <- range_bin[c(1,11)] %>% separate(bin, sep = "_", c('BIN_distance', 'BIN_CROSS.JTE'), remove = F)
    range_bin <- dcast(range_bin, BIN_distance~BIN_CROSS.JTE, value.var = "gene_N.distance_range.CROSS.JTE_range")
    range_bin <- range_bin[order(nrow(range_bin):1),]
    
    
    df2D0 <- df2D
    
    MedianAll <- median(df2D0$log2Ratio)
    df2D2 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(Median=median(log2Ratio))
    df2D2$value <- df2D2$Median - MedianAll
    
    df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_distance', 'BIN_CROSS.JTE'), remove = F)
    df2D3 <- dcast(df2D2, BIN_distance~BIN_CROSS.JTE, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_distance)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot[[paste0(cell,".",b)]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.log2", column_title = paste0(cell,".",b), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
  }
  
  
  DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
  names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
  CROSS_bin <- df_CROSS[c(1,11,12)]
  names(CROSS_bin) <- c("gene_symbol","distance","CROSS.JTE")
  DF <- merge(DGE,CROSS_bin)
  
  df2D <- DF
  binA1=binA2=binA3=binA4=round(nrow(df2D)/5)
  binA5=nrow(df2D)-(binA1+binA2+binA3+binA4)
  
  df2D <- df2D[order(df2D$distance,decreasing=F),]
  df2D$BIN_distance <- c(rep("BinD1",binA1),rep("BinD2",binA2),rep("BinD3",binA3),rep("BinD4",binA4),rep("BinD5",binA5))
  
  df2D <- df2D[order(df2D$CROSS.JTE,decreasing=F),]
  df2D$BIN_CROSS.JTE <- c(rep("binC1",binA1),rep("binC2",binA2),rep("binC3",binA3),rep("binC4",binA4),rep("binC5",binA5))
  
  
  df2D$bin <- paste0(df2D$BIN_distance,"_",df2D$BIN_CROSS.JTE)
  
  range_bin <- df2D %>%
    group_by(bin) %>%
    dplyr::summarise(gene_N_distance=length(distance),
                     min_distance=min(distance),
                     max_distance=max(distance),
                     gene_N_CROSS.JTE=length(CROSS.JTE),
                     min_CROSS.JTE=min(CROSS.JTE),
                     max_CROSS.JTE=max(CROSS.JTE))
  
  range_bin$distance_range <- paste(round(range_bin$min_distance,2),round(range_bin$max_distance,2),sep = " ~ ")
  range_bin$CROSS.JTE_range <- paste(round(range_bin$min_CROSS.JTE,2),round(range_bin$max_CROSS.JTE,2),sep = " ~ ")
  range_bin$gene_N <- range_bin$gene_N_distance
  range_bin$gene_N.distance_range.CROSS.JTE_range <- paste(range_bin$gene_N,range_bin$distance_range,range_bin$CROSS.JTE_range,sep = ";")
  range_bin <- range_bin[c(1,11)] %>% separate(bin, sep = "_", c('BIN_distance', 'BIN_CROSS.JTE'), remove = F)
  range_bin <- dcast(range_bin, BIN_distance~BIN_CROSS.JTE, value.var = "gene_N.distance_range.CROSS.JTE_range")
  range_bin <- range_bin[order(nrow(range_bin):1),]
  
  
  df2D0 <- df2D
  
  MedianAll <- median(df2D0$log2Ratio)
  df2D2 <- df2D0 %>%
    group_by(bin) %>%
    dplyr::summarise(Median=median(log2Ratio))
  df2D2$value <- df2D2$Median# - MedianAll
  
  df2D2 <- df2D2 %>% separate(bin, sep = "_", c('BIN_distance', 'BIN_CROSS.JTE'), remove = F)
  df2D3 <- dcast(df2D2, BIN_distance~BIN_CROSS.JTE, value.var = "value")
  matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_distance)
  matrix2D <- matrix2D[order(nrow(matrix2D):1),]
  matrix2D <- as.matrix(matrix2D)
  
  hm_plot_all[[paste0(cell)]] <- 
    grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.log2", column_title = paste0(cell), 
                               row_order = rownames(matrix2D), 
                               column_order = colnames(matrix2D),
                               col = col_fun)))
  
}

library(cowplot)
# plot_grid(plotlist=scatterplot,ncol = 5)
plot_grid(plotlist=hm_plot,ncol = 5)
plot_grid(plotlist=hm_plot[c(16:20)],ncol = 5)

plot_grid(plotlist=hm_plot_all,ncol = 4)
plot_grid(plotlist=hm_plot_all[c(4)],ncol = 1)


# write.csv(matrix2D, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig4e.csv"),row.names = T)
