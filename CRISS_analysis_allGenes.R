workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

library(tidyverse)
library(RColorBrewer)

DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))
DNG_sd <- DNGlist[["sd"]][c(4,9)]
names(DNG_sd)[c(1)] <- "gene_symbol"

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CRISS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CRISS_table/List.CRISS.RData"))

List.CRISSa <- List.CRISS

List.CRISSa <- lapply(List.CRISSa, function(x) { x <- merge(DNG_sd,x,by.x = "gene_symbol",by.y = "gene_symbol"); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })
# 
# DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
# names(DRas_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, nb_gene_symbol %in% DRas_clean_genes$gene_symbol); x })


ggplot() +
  stat_ecdf(data=List.CRISSa[[c(1)]], aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSa[[c(2)]], aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSa[[c(3)]], aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSa[[c(4)]], aes(x=CRISS_JTE.v.DMSO, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
                              brewer.pal(n = 12,name = 'Paired')),
                     labels = c(paste0(names(List.CRISSa)[c(1)],", ",nrow(List.CRISSa[[c(1)]]),", ",round(median(List.CRISSa[[c(1)]][[c(6)]]),2)),
                                paste0(names(List.CRISSa)[c(2)],", ",nrow(List.CRISSa[[c(2)]]),", ",round(median(List.CRISSa[[c(2)]][[c(6)]]),2)),
                                paste0(names(List.CRISSa)[c(3)],", ",nrow(List.CRISSa[[c(3)]]),", ",round(median(List.CRISSa[[c(3)]][[c(6)]]),2)),
                                paste0(names(List.CRISSa)[c(4)],", ",nrow(List.CRISSa[[c(4)]]),", ",round(median(List.CRISSa[[c(4)]][[c(6)]]),2)))) +
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
  coord_cartesian(xlim = c(0, 3)) +
  labs(title = paste("CRISS, all UR-isolated genes",sep = ""),
       x = "CRISS", y = "Cumulative fraction", color = "")





List.DGE_QSF2 <- lapply(List.DGE_QSF, function(x) { x <- x[c(1,4,5)]; x })
DGE_QSF <- Reduce(function(x, y) merge(x, y), List.DGE_QSF2)
DGE_QSF2 <- subset(DGE_QSF, DGE_QSF[[2]]>=10 & DGE_QSF[[3]]>=10 & DGE_QSF[[4]]>=10 & DGE_QSF[[5]]>=10 & DGE_QSF[[6]]>=10 & DGE_QSF[[7]]>=10 & DGE_QSF[[8]]>=10 & DGE_QSF[[9]]>=10)


List.CRISSb <- lapply(List.CRISSa, function(x) { x <- subset(x,gene_symbol %in% DGE_QSF2$gene_symbol); x })

ggplot() +
  stat_ecdf(data=List.CRISSb[[c(1)]], aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSb[[c(2)]], aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSb[[c(3)]], aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=List.CRISSb[[c(4)]], aes(x=CRISS_JTE.v.DMSO, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
    brewer.pal(n = 12,name = 'Paired')),
    labels = c(paste0(names(List.CRISSb)[c(1)],", ",nrow(List.CRISSb[[c(1)]]),", ",round(median(List.CRISSb[[c(1)]][[c(6)]]),2)),
               paste0(names(List.CRISSb)[c(2)],", ",nrow(List.CRISSb[[c(2)]]),", ",round(median(List.CRISSb[[c(2)]][[c(6)]]),2)),
               paste0(names(List.CRISSb)[c(3)],", ",nrow(List.CRISSb[[c(3)]]),", ",round(median(List.CRISSb[[c(3)]][[c(6)]]),2)),
               paste0(names(List.CRISSb)[c(4)],", ",nrow(List.CRISSb[[c(4)]]),", ",round(median(List.CRISSb[[c(4)]][[c(6)]]),2)))) +
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
  coord_cartesian(xlim = c(0, 3)) +
  labs(title = paste("CRISS, commonly expressed UR-isolated genes (RPM >=10 in all QSF samples)",sep = ""),
       x = "CRISS", y = "Cumulative fraction", color = "")




#### GOI_QSF_DGE vs CRISS


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
DNG_sd <- DNGlist[["sd"]][c(4,9)]
names(DNG_sd)[c(1)] <- "gene_symbol"

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CRISS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CRISS_table/List.CRISS.RData"))

List.CRISSa <- List.CRISS

List.CRISSa <- lapply(List.CRISSa, function(x) { x <- merge(DNG_sd,x,by.x = "gene_symbol",by.y = "gene_symbol"); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })
# 
# DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
# names(DRas_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, nb_gene_symbol %in% DRas_clean_genes$gene_symbol); x })

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
scatterplot <- list()

CRISS_plotTbl <-list()
CRISS_plotTbl_out <-list()
CRISS_plotTbl_out2 <-list()
CRISS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
KS_test <- list()

for (cell in cells) {

  # cell="HepG2"
  
  ldf <- List.DGE_QSF[[cell]]
  QSF_UP <- subset(ldf,ldf[[10]]=="UP")
  QSF_NC <- subset(ldf,ldf[[10]]=="NO")
  QSF_DN <- subset(ldf,ldf[[10]]=="DN")
  
  
  CRISS <- List.CRISSa[[cell]]

  CDFplot[[cell]] <- ggplot() +
    stat_ecdf(data=subset(CRISS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(CRISS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(CRISS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(cell, ", QSF.DGE vs delta.CRISS",
                        "\nNCvUP, ", formatC(ks.test(subset(CRISS, gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO, subset(CRISS, gene_symbol %in% QSF_UP$gene_symbol)$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2),
                        "\nNCvDN, ", formatC(ks.test(subset(CRISS, gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO, subset(CRISS, gene_symbol %in% QSF_DN$gene_symbol)$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2)),
         x = "CRISS", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(CRISS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(CRISS, gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO),2)), 
                                  paste0("UP genes, ",nrow(subset(CRISS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(CRISS, gene_symbol %in% QSF_UP$gene_symbol)$CRISS_JTE.v.DMSO),2)),
                                  paste0("DN genes, ",nrow(subset(CRISS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(CRISS, gene_symbol %in% QSF_DN$gene_symbol)$CRISS_JTE.v.DMSO),2))
                                  
    ), 
    values=c("grey","red","blue")) + 
    coord_cartesian(xlim = c(-1, 3))
    
  
  DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
  names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
  CRISS <- List.CRISSa[[cell]]
  DF <- merge(DGE,CRISS)
  
  cor_all <- cor.test(DF$log2Ratio,DF$CRISS_JTE.v.DMSO)
  scatterplot[[cell]] <- ggplot(data=DF, aes(x=log2Ratio, y=CRISS_JTE.v.DMSO, label=gene_symbol)) +
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
         x = "log2Ratio", y = "CRISS_JTE.v.DMSO")
  
  
  
  ######### 5 bins by GOI expression level in DMSO 
  
  DGE <- List.DGE_QSF[[cell]][c(1,4,6,10)]
  names(DGE) <- c("gene_symbol","RPM_DMSO","log2Ratio","regu.adj")
  CRISS <- List.CRISSa[[cell]]
  df_CRISS <- merge(DGE,CRISS)
  
  df_CRISS <- df_CRISS[order(df_CRISS$RPM_DMSO,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CRISS)/5)
  bin5=nrow(df_CRISS)-(bin1+bin2+bin3+bin4)
  df_CRISS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CRISS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CRISS, bin == "bin1")),
                                        nrow(subset(df_CRISS, bin == "bin2")),
                                        nrow(subset(df_CRISS, bin == "bin3")),
                                        nrow(subset(df_CRISS, bin == "bin4")),
                                        nrow(subset(df_CRISS, bin == "bin5"))),
                             "RPM_DMSO_range"=c(paste(range(subset(df_CRISS, bin == "bin1")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin2")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin3")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin4")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin5")$RPM_DMSO),collapse = "~")))
  names(CRISS_bin_tbl2) <- paste(names(CRISS_bin_tbl2),cell,sep = "_")
  CRISS_bin_tbl <- cbind(CRISS_bin_tbl,CRISS_bin_tbl2)
  
  for (b in c("bin1","bin2","bin3","bin4","bin5")){
    
    # b="bin1"
    
    ldf <- subset(df_CRISS,bin==b)
    
    QSF_UP <- subset(ldf,regu.adj=="UP")
    QSF_NC <- subset(ldf,regu.adj=="NO")
    QSF_DN <- subset(ldf,regu.adj=="DN")
    
    
    CDFplot[[paste0(cell,"_RPM_DMSO_",b)]] <- ggplot() +
      stat_ecdf(data=QSF_NC, aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
      stat_ecdf(data=QSF_UP, aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
      stat_ecdf(data=QSF_DN, aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
      theme_bw() +
      theme(plot.title = element_text(size = 10)) +
      labs(title = paste0(cell,"_RPM_DMSO_",b, ", QSF.DGE vs CRISS",
                          "\nNCvUP, ", formatC(ks.test(QSF_NC$CRISS_JTE.v.DMSO, QSF_UP$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2),
                          "\nNCvDN, ", formatC(ks.test(QSF_NC$CRISS_JTE.v.DMSO, QSF_DN$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2)),
           x = "CRISS", y = "Cumulative fraction", color = "") +
      scale_color_manual(labels = c(paste0("NC genes, ",nrow(QSF_NC),", ",round(median(QSF_NC$CRISS_JTE.v.DMSO),2)), 
                                    paste0("UP genes, ",nrow(QSF_UP),", ",round(median(QSF_UP$CRISS_JTE.v.DMSO),2)),
                                    paste0("DN genes, ",nrow(QSF_DN),", ",round(median(QSF_DN$CRISS_JTE.v.DMSO),2))
                                    
      ), 
      values=c("grey","red","blue")) + 
      coord_cartesian(xlim = c(-1, 3))
    
    
  }
}

library(cowplot)
plot_grid(plotlist=CDFplot,ncol = 6)
plot_grid(plotlist=CDFplot[c(19)])
plot_grid(plotlist=CDFplot[c(20:24)],ncol = 3)

# plot_grid(plotlist=scatterplot,ncol = 4)



#### NNG_QSF_DGE vs CRISS


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
DNG_sd <- DNGlist[["sd"]][c(4,9)]
names(DNG_sd)[c(1)] <- "gene_symbol"

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CRISS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CRISS_table/List.CRISS.RData"))

List.CRISSa <- List.CRISS

List.CRISSa <- lapply(List.CRISSa, function(x) { x <- merge(DNG_sd,x,by.x = "gene_symbol",by.y = "gene_symbol"); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })
# 
# DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
# names(DRas_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, nb_gene_symbol %in% DRas_clean_genes$gene_symbol); x })

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
scatterplot <- list()

CRISS_plotTbl <-list()
CRISS_plotTbl_out <-list()
CRISS_plotTbl_out2 <-list()
CRISS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))
KS_test <- list()

for (cell in cells) {
  
  # cell="HepG2"
  
  ldf <- List.DGE_QSF[[cell]]
  QSF_UP <- subset(ldf,ldf[[10]]=="UP")
  QSF_NC <- subset(ldf,ldf[[10]]=="NO")
  QSF_DN <- subset(ldf,ldf[[10]]=="DN")
  
  
  CRISS <- List.CRISSa[[cell]]
  
  CDFplot[[cell]] <- ggplot() +
    stat_ecdf(data=subset(CRISS, nb_gene_symbol %in% QSF_NC$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(CRISS, nb_gene_symbol %in% QSF_UP$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(CRISS, nb_gene_symbol %in% QSF_DN$gene_symbol), aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(cell, ", QSF.DGE vs delta.CRISS",
                        "\nNCvUP, ", formatC(ks.test(subset(CRISS, nb_gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO, subset(CRISS, nb_gene_symbol %in% QSF_UP$gene_symbol)$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2),
                        "\nNCvDN, ", formatC(ks.test(subset(CRISS, nb_gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO, subset(CRISS, nb_gene_symbol %in% QSF_DN$gene_symbol)$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2)),
         x = "CRISS", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(CRISS, nb_gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(CRISS, nb_gene_symbol %in% QSF_NC$gene_symbol)$CRISS_JTE.v.DMSO),2)), 
                                  paste0("UP genes, ",nrow(subset(CRISS, nb_gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(CRISS, nb_gene_symbol %in% QSF_UP$gene_symbol)$CRISS_JTE.v.DMSO),2)),
                                  paste0("DN genes, ",nrow(subset(CRISS, nb_gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(CRISS, nb_gene_symbol %in% QSF_DN$gene_symbol)$CRISS_JTE.v.DMSO),2))
                                  
    ), 
    values=c("grey","red","blue")) + 
    coord_cartesian(xlim = c(-1, 3))
  
  # write.csv(merge(CRISS[c(2,7)],ldf[c(1,10)],by.x = "nb_gene_symbol",by.y = "gene_symbol"), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig4k_2.csv"),row.names = T)
  
  
  DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
  names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
  CRISS <- List.CRISSa[[cell]]
  DF <- merge(DGE,CRISS,by.x = "gene_symbol",by.y = "nb_gene_symbol")
  
  cor_all <- cor.test(DF$log2Ratio,DF$CRISS_JTE.v.DMSO)
  scatterplot[[cell]] <- ggplot(data=DF, aes(x=log2Ratio, y=CRISS_JTE.v.DMSO, label=gene_symbol)) +
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
         x = "log2Ratio", y = "CRISS_JTE.v.DMSO")
  
  
  
  
  ######### 5 bins by GOI expression level in DMSO 
  
  DGE <- List.DGE_QSF[[cell]][c(1,4,6,10)]
  names(DGE) <- c("gene_symbol","RPM_DMSO","log2Ratio","regu.adj")
  CRISS <- List.CRISSa[[cell]]
  df_CRISS <- merge(DGE,CRISS,by.x = "gene_symbol",by.y = "nb_gene_symbol")
  
  df_CRISS <- df_CRISS[order(df_CRISS$RPM_DMSO,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CRISS)/5)
  bin5=nrow(df_CRISS)-(bin1+bin2+bin3+bin4)
  df_CRISS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CRISS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CRISS, bin == "bin1")),
                                        nrow(subset(df_CRISS, bin == "bin2")),
                                        nrow(subset(df_CRISS, bin == "bin3")),
                                        nrow(subset(df_CRISS, bin == "bin4")),
                                        nrow(subset(df_CRISS, bin == "bin5"))),
                             "RPM_DMSO_range"=c(paste(range(subset(df_CRISS, bin == "bin1")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin2")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin3")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin4")$RPM_DMSO),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin5")$RPM_DMSO),collapse = "~")))
  names(CRISS_bin_tbl2) <- paste(names(CRISS_bin_tbl2),cell,sep = "_")
  CRISS_bin_tbl <- cbind(CRISS_bin_tbl,CRISS_bin_tbl2)
  
  for (b in c("bin1","bin2","bin3","bin4","bin5")){
    
    # b="bin1"
    
    ldf <- subset(df_CRISS,bin==b)
    
    QSF_UP <- subset(ldf,regu.adj=="UP")
    QSF_NC <- subset(ldf,regu.adj=="NO")
    QSF_DN <- subset(ldf,regu.adj=="DN")
    
    
    CDFplot[[paste0(cell,"_RPM_DMSO_",b)]] <- ggplot() +
      stat_ecdf(data=QSF_NC, aes(x=CRISS_JTE.v.DMSO, color = "a"),geom = "step", size=1) +
      stat_ecdf(data=QSF_UP, aes(x=CRISS_JTE.v.DMSO, color = "b"),geom = "step", size=1) +
      stat_ecdf(data=QSF_DN, aes(x=CRISS_JTE.v.DMSO, color = "c"),geom = "step", size=1) +
      theme_bw() +
      theme(plot.title = element_text(size = 10)) +
      labs(title = paste0(cell,"_RPM_DMSO_",b, ", QSF.DGE vs CRISS",
                          "\nNCvUP, ", formatC(ks.test(QSF_NC$CRISS_JTE.v.DMSO, QSF_UP$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2),
                          "\nNCvDN, ", formatC(ks.test(QSF_NC$CRISS_JTE.v.DMSO, QSF_DN$CRISS_JTE.v.DMSO)$p.value, format = "e", digits = 2)),
           x = "CRISS", y = "Cumulative fraction", color = "") +
      scale_color_manual(labels = c(paste0("NC genes, ",nrow(QSF_NC),", ",round(median(QSF_NC$CRISS_JTE.v.DMSO),2)), 
                                    paste0("UP genes, ",nrow(QSF_UP),", ",round(median(QSF_UP$CRISS_JTE.v.DMSO),2)),
                                    paste0("DN genes, ",nrow(QSF_DN),", ",round(median(QSF_DN$CRISS_JTE.v.DMSO),2))
                                    
      ), 
      values=c("grey","red","blue")) + 
      coord_cartesian(xlim = c(-1, 3))
    
    # write.csv(ldf[c(1,4,10)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig4k_1.csv"),row.names = T)
    
  }
}

library(cowplot)
plot_grid(plotlist=CDFplot,ncol = 6)
plot_grid(plotlist=CDFplot[c(19)])
plot_grid(plotlist=CDFplot[c(20)])
plot_grid(plotlist=CDFplot[c(20:24)],ncol = 3)

plot_grid(plotlist=scatterplot,ncol = 4)




#### CRISS vs distance

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
DNG_sd <- DNGlist[["sd"]][c(4,9,11)]
names(DNG_sd)[c(1)] <- "gene_symbol"

plot <- list()
CRISS_plotTbl <-list()
CRISS_plotTbl_out <-list()
CRISS_plotTbl_out2 <-list()

CRISS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

List.CRISS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CRISS_table/List.CRISS.RData"))

List.CRISSa <- List.CRISS

List.CRISSa <- lapply(List.CRISSa, function(x) { x <- merge(DNG_sd,x,by.x = "gene_symbol",by.y = "gene_symbol"); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })
# 
# DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
# names(DRas_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, nb_gene_symbol %in% DRas_clean_genes$gene_symbol); x })



cells <- c("UP", "U937", "HeLa", "HepG2")

KS_test <- list()

for (cell in cells) {

  # cell="HepG2"

  df_CRISS <- List.CRISSa[[cell]]
  names(df_CRISS)[c(8)] <- c("CRISS")


  df_CRISS <- df_CRISS[order(df_CRISS$distance,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CRISS)/5)
  bin5=nrow(df_CRISS)-(bin1+bin2+bin3+bin4)
  df_CRISS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))

  CRISS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CRISS, bin == "bin1")),
                                      nrow(subset(df_CRISS, bin == "bin2")),
                                      nrow(subset(df_CRISS, bin == "bin3")),
                                      nrow(subset(df_CRISS, bin == "bin4")),
                                      nrow(subset(df_CRISS, bin == "bin5"))),
                           "distance_range"=c(paste(range(subset(df_CRISS, bin == "bin1")$distance),collapse = "~"),
                                              paste(range(subset(df_CRISS, bin == "bin2")$distance),collapse = "~"),
                                              paste(range(subset(df_CRISS, bin == "bin3")$distance),collapse = "~"),
                                              paste(range(subset(df_CRISS, bin == "bin4")$distance),collapse = "~"),
                                              paste(range(subset(df_CRISS, bin == "bin5")$distance),collapse = "~")))
  names(CRISS_bin_tbl2) <- paste(names(CRISS_bin_tbl2),cell,sep = "_")
  CRISS_bin_tbl <- cbind(CRISS_bin_tbl,CRISS_bin_tbl2)

  df_CRISS3 <- data_summary(df_CRISS, varname="CRISS",
                          groupnames=c("bin"))
  df_CRISS3$sampleName <- rep(cell, nrow(df_CRISS3))
  CRISS_plotTbl[[cell]] <- df_CRISS3

  df_CRISS3_out <- df_CRISS3
  df_CRISS3_out$sampleName <- rep(paste(cell,sep = "_"), nrow(df_CRISS3_out))
  CRISS_plotTbl_out[[cell]] <- df_CRISS3_out

  KS_test[[cell]] <- ks.test(subset(df_CRISS, bin == "bin1")$CRISS,subset(df_CRISS, bin == "bin5")$CRISS)



}


CRISS_plotTbl_JTE <- do.call(rbind,CRISS_plotTbl)
CRISS_plotTbl_out2 <- do.call(rbind,CRISS_plotTbl_out)
# write.csv(CRISS_plotTbl_out2, "./CRISS_plotTbl.csv", col.names = F, row.names = T)

library(RColorBrewer)
ggplot(CRISS_plotTbl_JTE, aes(x=bin, y=CRISS, group=sampleName, color=sampleName)) +
  geom_line(size=1) +
  geom_point(size=3) +
  # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
  geom_errorbar(aes(ymin=CRISS-se, ymax=CRISS+se), width=0.5, size=1,
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
       x = "distance size bin", y = "Mean of CRISS") +
  scale_color_manual(values=c(#"gray80","gray50","gray0",
    brewer.pal(n = 12,name = 'Paired')[c(1,2,3,4)]))











#### DGE vs CRISS by distance

library(tidyverse)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

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

plot <- list()
CRISS_plotTbl <-list()
CRISS_plotTbl_out <-list()
CRISS_plotTbl_out2 <-list()

CRISS_bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))


DNGlist <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/DNG/DNGlist_20230317.RData"))
DNG_sd <- DNGlist[["sd"]][c(4,9,11)]
names(DNG_sd)[c(1)] <- "gene_symbol"

List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))[c(2,4,5,6)]
names(List.DGE_QSF) <- c("HeLa", "HepG2", "U937", "UP")

List.CRISS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_CRISS_table/List.CRISS.RData"))

List.CRISSa <- List.CRISS

List.CRISSa <- lapply(List.CRISSa, function(x) { x <- merge(DNG_sd,x,by.x = "gene_symbol",by.y = "gene_symbol"); x })

# URs_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/URs.clean.genes_20230330.tbl"), header = F, quote = "")
# names(URs_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, gene_symbol %in% URs_clean_genes$gene_symbol); x })
# 
# DRas_clean_genes <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/_clean_genes/DRas.clean.genes_20230330.tbl"), header = F, quote = "")
# names(DRas_clean_genes) <- "gene_symbol"
# List.CRISSa <- lapply(List.CRISSa, function(x) { x <- subset(x, nb_gene_symbol %in% DRas_clean_genes$gene_symbol); x })


cells <- c("UP", "U937", "HeLa", "HepG2")

KS_test <- list()
scatterplot_GOI <- list()
scatterplot_NNG <- list()
hm_plot_GOI <- list()
hm_plot_NNG <- list()

for (cell in cells) {
  
  # cell="UP"
  
  df_CRISS <- List.CRISSa[[cell]][c(1,2,3,8)]
  names(df_CRISS)[c(4)] <- "CRISS"
  
  
  df_CRISS <- df_CRISS[order(df_CRISS$distance,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(df_CRISS)/5)
  bin5=nrow(df_CRISS)-(bin1+bin2+bin3+bin4)
  df_CRISS$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  CRISS_bin_tbl2<-data.frame("gene_N"=c(nrow(subset(df_CRISS, bin == "bin1")),
                                        nrow(subset(df_CRISS, bin == "bin2")),
                                        nrow(subset(df_CRISS, bin == "bin3")),
                                        nrow(subset(df_CRISS, bin == "bin4")),
                                        nrow(subset(df_CRISS, bin == "bin5"))),
                             "distance_range"=c(paste(range(subset(df_CRISS, bin == "bin1")$distance),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin2")$distance),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin3")$distance),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin4")$distance),collapse = "~"),
                                                paste(range(subset(df_CRISS, bin == "bin5")$distance),collapse = "~")))
  names(CRISS_bin_tbl2) <- paste(names(CRISS_bin_tbl2),cell,sep = "_")
  CRISS_bin_tbl <- cbind(CRISS_bin_tbl,CRISS_bin_tbl2)
  
  df_CRISS3 <- data_summary(df_CRISS, varname="CRISS", 
                            groupnames=c("bin"))
  df_CRISS3$sampleName <- rep(cell, nrow(df_CRISS3))
  CRISS_plotTbl[[cell]] <- df_CRISS3
  
  df_CRISS3_out <- df_CRISS3
  df_CRISS3_out$sampleName <- rep(paste(cell,sep = "_"), nrow(df_CRISS3_out))
  CRISS_plotTbl_out[[cell]] <- df_CRISS3_out
  
  KS_test[[cell]] <- ks.test(subset(df_CRISS, bin == "bin1")$CRISS,subset(df_CRISS, bin == "bin5")$CRISS)
  
  
  for (b in c("bin1","bin2","bin3","bin4","bin5")){
    
    # b="bin1"
    
    DGE <- List.DGE_QSF[[cell]][c(1,6,10)]
    names(DGE) <- c("gene_symbol","log2Ratio","regu.adj")
    
    
    CRISS_bin <- subset(df_CRISS,bin==b)[c(1,4)]
    names(CRISS_bin) <- c("gene_symbol","CRISS_JTE.v.DMSO")
    DF <- merge(DGE,CRISS_bin)
    
    cor_all <- cor.test(DF$log2Ratio,DF$CRISS_JTE.v.DMSO)
    scatterplot_GOI[[paste0(cell,"_",b)]] <- ggplot(data=DF, aes(x=log2Ratio, y=CRISS_JTE.v.DMSO, label=gene_symbol)) +
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
           x = "log2Ratio", y = "CRISS_JTE.v.DMSO")
    
    
    
    df2D <- DF
    binA1=binA2=binA3=binA4=round(nrow(df2D)/5)
    binA5=nrow(df2D)-(binA1+binA2+binA3+binA4)
    
    df2D <- df2D[order(df2D$log2Ratio,decreasing=F),]
    df2D$BIN_log2Ratio <- c(rep("BinL1",binA1),rep("BinL2",binA2),rep("BinL3",binA3),rep("BinL4",binA4),rep("BinL5",binA5))
    
    df2D <- df2D[order(df2D$CRISS_JTE.v.DMSO,decreasing=F),]
    df2D$BIN_CRISS_JTE.v.DMSO <- c(rep("binC1",binA1),rep("binC2",binA2),rep("binC3",binA3),rep("binC4",binA4),rep("binC5",binA5))
    
    
    df2D$bin <- paste0(df2D$BIN_log2Ratio,"_",df2D$BIN_CRISS_JTE.v.DMSO)
    
    range_bin <- df2D %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_log2Ratio=length(log2Ratio),
                       min_log2Ratio=min(log2Ratio),
                       max_log2Ratio=max(log2Ratio),
                       gene_N_CRISS_JTE.v.DMSO=length(CRISS_JTE.v.DMSO),
                       min_CRISS_JTE.v.DMSO=min(CRISS_JTE.v.DMSO),
                       max_CRISS_JTE.v.DMSO=max(CRISS_JTE.v.DMSO))
    
    range_bin$log2Ratio_range <- paste(round(range_bin$min_log2Ratio,2),round(range_bin$max_log2Ratio,2),sep = " ~ ")
    range_bin$CRISS_JTE.v.DMSO_range <- paste(round(range_bin$min_CRISS_JTE.v.DMSO,2),round(range_bin$max_CRISS_JTE.v.DMSO,2),sep = " ~ ")
    range_bin$gene_N <- range_bin$gene_N_log2Ratio
    range_bin$gene_N.log2Ratio_range.CRISS_JTE.v.DMSO_range <- paste(range_bin$gene_N,range_bin$log2Ratio_range,range_bin$CRISS_JTE.v.DMSO_range,sep = ";")
    range_bin <- range_bin[c(1,11)] %>% separate(bin, sep = "_", c('BIN_log2Ratio', 'BIN_CRISS_JTE.v.DMSO'), remove = F)
    range_bin <- dcast(range_bin, BIN_log2Ratio~BIN_CRISS_JTE.v.DMSO, value.var = "gene_N.log2Ratio_range.CRISS_JTE.v.DMSO_range")
    range_bin <- range_bin[order(nrow(range_bin):1),]
    
    
    df2D0 <- df2D
    
    range_bin0 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_log2Ratio=length(log2Ratio),
                       min_log2Ratio=min(log2Ratio),
                       max_log2Ratio=max(log2Ratio),
                       gene_N_CRISS_JTE.v.DMSO=length(CRISS_JTE.v.DMSO),
                       min_CRISS_JTE.v.DMSO=min(CRISS_JTE.v.DMSO),
                       max_CRISS_JTE.v.DMSO=max(CRISS_JTE.v.DMSO))
    range_bin0$gene_N <- range_bin0$gene_N_log2Ratio
    
    
    range_bin0$NormValue <- range_bin0$gene_N/sum(range_bin0$gene_N)
    MedianAll <- median(range_bin0$NormValue)
    range_bin0$value <- log2(range_bin0$NormValue/MedianAll)
    
    range_bin0 <- range_bin0 %>% separate(bin, sep = "_", c('BIN_log2Ratio', 'BIN_CRISS_JTE.v.DMSO'), remove = F)
    df2D3 <- dcast(range_bin0, BIN_log2Ratio~BIN_CRISS_JTE.v.DMSO, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_log2Ratio)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot_GOI[[paste0(cell,".",b)]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.log2", column_title = paste0(cell,".",b), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
    
    
    
    
    
    
    CRISS_bin <- subset(df_CRISS,bin==b)[c(2,4)]
    names(CRISS_bin) <- c("gene_symbol","CRISS_JTE.v.DMSO")
    DF <- merge(DGE,CRISS_bin)
    
    cor_all <- cor.test(DF$log2Ratio,DF$CRISS_JTE.v.DMSO)
    scatterplot_NNG[[paste0(cell,"_",b)]] <- ggplot(data=DF, aes(x=log2Ratio, y=CRISS_JTE.v.DMSO, label=gene_symbol)) +
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
           x = "log2Ratio", y = "CRISS_JTE.v.DMSO")
    
    
    
    df2D <- DF
    binA1=binA2=binA3=binA4=round(nrow(df2D)/5)
    binA5=nrow(df2D)-(binA1+binA2+binA3+binA4)
    
    df2D <- df2D[order(df2D$log2Ratio,decreasing=F),]
    df2D$BIN_log2Ratio <- c(rep("BinL1",binA1),rep("BinL2",binA2),rep("BinL3",binA3),rep("BinL4",binA4),rep("BinL5",binA5))
    
    df2D <- df2D[order(df2D$CRISS_JTE.v.DMSO,decreasing=F),]
    df2D$BIN_CRISS_JTE.v.DMSO <- c(rep("binC1",binA1),rep("binC2",binA2),rep("binC3",binA3),rep("binC4",binA4),rep("binC5",binA5))
    
    
    df2D$bin <- paste0(df2D$BIN_log2Ratio,"_",df2D$BIN_CRISS_JTE.v.DMSO)
    
    range_bin <- df2D %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_log2Ratio=length(log2Ratio),
                       min_log2Ratio=min(log2Ratio),
                       max_log2Ratio=max(log2Ratio),
                       gene_N_CRISS_JTE.v.DMSO=length(CRISS_JTE.v.DMSO),
                       min_CRISS_JTE.v.DMSO=min(CRISS_JTE.v.DMSO),
                       max_CRISS_JTE.v.DMSO=max(CRISS_JTE.v.DMSO))
    
    range_bin$log2Ratio_range <- paste(round(range_bin$min_log2Ratio,2),round(range_bin$max_log2Ratio,2),sep = " ~ ")
    range_bin$CRISS_JTE.v.DMSO_range <- paste(round(range_bin$min_CRISS_JTE.v.DMSO,2),round(range_bin$max_CRISS_JTE.v.DMSO,2),sep = " ~ ")
    range_bin$gene_N <- range_bin$gene_N_log2Ratio
    range_bin$gene_N.log2Ratio_range.CRISS_JTE.v.DMSO_range <- paste(range_bin$gene_N,range_bin$log2Ratio_range,range_bin$CRISS_JTE.v.DMSO_range,sep = ";")
    range_bin <- range_bin[c(1,11)] %>% separate(bin, sep = "_", c('BIN_log2Ratio', 'BIN_CRISS_JTE.v.DMSO'), remove = F)
    range_bin <- dcast(range_bin, BIN_log2Ratio~BIN_CRISS_JTE.v.DMSO, value.var = "gene_N.log2Ratio_range.CRISS_JTE.v.DMSO_range")
    range_bin <- range_bin[order(nrow(range_bin):1),]
    
    
    df2D0 <- df2D
    
    range_bin0 <- df2D0 %>%
      group_by(bin) %>%
      dplyr::summarise(gene_N_log2Ratio=length(log2Ratio),
                       min_log2Ratio=min(log2Ratio),
                       max_log2Ratio=max(log2Ratio),
                       gene_N_CRISS_JTE.v.DMSO=length(CRISS_JTE.v.DMSO),
                       min_CRISS_JTE.v.DMSO=min(CRISS_JTE.v.DMSO),
                       max_CRISS_JTE.v.DMSO=max(CRISS_JTE.v.DMSO))
    range_bin0$gene_N <- range_bin0$gene_N_log2Ratio
    
    
    range_bin0$NormValue <- range_bin0$gene_N/sum(range_bin0$gene_N)
    MedianAll <- median(range_bin0$NormValue)
    range_bin0$value <- log2(range_bin0$NormValue/MedianAll)
    
    range_bin0 <- range_bin0 %>% separate(bin, sep = "_", c('BIN_log2Ratio', 'BIN_CRISS_JTE.v.DMSO'), remove = F)
    df2D3 <- dcast(range_bin0, BIN_log2Ratio~BIN_CRISS_JTE.v.DMSO, value.var = "value")
    matrix2D <- data.frame(df2D3[-c(1)],row.names=df2D3$BIN_log2Ratio)
    matrix2D <- matrix2D[order(nrow(matrix2D):1),]
    matrix2D <- as.matrix(matrix2D)
    
    hm_plot_NNG[[paste0(cell,".",b)]] <- 
      grid.grabExpr(draw(Heatmap(matrix2D, name = "Norm.log2", column_title = paste0(cell,".",b), 
                                 row_order = rownames(matrix2D), 
                                 column_order = colnames(matrix2D),
                                 col = col_fun)))
    
    
    
  }

}

library(cowplot)
plot_grid(plotlist=scatterplot_GOI,ncol = 5)
plot_grid(plotlist=scatterplot_NNG,ncol = 5)

plot_grid(plotlist=hm_plot_GOI,ncol = 5)
plot_grid(plotlist=hm_plot_NNG,ncol = 5)
