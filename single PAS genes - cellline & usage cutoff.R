
#####################################################



#### using genes with only 1 PAS

# The log2(RPM of FT sample/RPM of 4sU sample) value is called Stability Score (SS). Two biological replicates were averaged.

library(reshape2)
library(tidyverse)
library(ggrepel)

plot <- list()
for (cell in c("HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h","U937_1_7h","U937_PMA_1_7h")) {
  
  # cell <- "HepG2_10_7h"
  
  if (cell == "HeLa_1_9h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/HeLa/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(4)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  if (cell == "HeLa_10_9h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/HeLa/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(5)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  if (cell == "HepG2_1_7h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/HepG2/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(10)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  if (cell == "HepG2_10_7h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/HepG2/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(11)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  if (cell == "U937_1_7h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/U937_6h/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(15)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  if (cell == "U937_PMA_1_7h") {
    PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/U937_PMA_6h/pas_17.4.csv"), header = T)
    PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
    
    filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(16)]
    filenames
    DGE<-read.table(filenames, header = T, sep = ",")[c(1,6,10)]
    names(DGE)[c(2)] <- "log2"
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  PAS_N1 <- subset(PAS_aaa0, region=="3UTR" & sumCount >= 95) %>% count(gene_symbol)
  PAS_aaa1 <- subset(subset(PAS_aaa0, region=="3UTR" & sumCount >= 95), gene_symbol %in% subset(PAS_N1, n=1)$gene_symbol)
  
  
  PAS_N2 <- subset(PAS_aaa0, region=="3UTR" & sumCount >= 10) %>% count(gene_symbol)
  PAS_aaa2 <- subset(subset(PAS_aaa0, region=="3UTR" & sumCount >= 10), gene_symbol %in% subset(PAS_N2, n>=2)$gene_symbol)
  
  
  library(RColorBrewer)
  
  plot[[cell]] <- ggplot() +
    stat_ecdf(data=DGE, aes(x=log2, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol), aes(x=log2, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol), aes(x=log2, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    labs(
      title = paste0("DGE, ",cell,", PAS number is based on 3'UTR",
                     "\nAllvSPA, ", formatC(ks.test(DGE$log2, subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$log2)$p.value, format = "e", digits = 2),
                     "\nAllvAPA, ", formatC(ks.test(DGE$log2, subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$log2)$p.value, format = "e", digits = 2),
                     "\nAPAvSPA, ", formatC(ks.test(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$log2, subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$log2)$p.value, format = "e", digits = 2)
      ),
      x = "log2FoldChange", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("all genes, ",nrow(DGE),", ",round(median(DGE$log2),2)), 
                                  paste0("SPA genes (>=95% usage), ",nrow(subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)),", ",round(median(subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$log2),2)), 
                                  paste0(">=2PAS genes (>=10% usage), ",nrow(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)),", ",round(median(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$log2),2))),
                       values=c("grey",brewer.pal(n = 8,name = 'Dark2'))) +
    scale_x_continuous(limits = c(-2, 1))
  
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)

DGE$Type_by_PAS_number <- "other"
DGE$Type_by_PAS_number[DGE$gene_symbol %in% PAS_aaa1$gene_symbol] <- "SPA"
DGE$Type_by_PAS_number[DGE$gene_symbol %in% PAS_aaa2$gene_symbol] <- ">=2PAS"
# write.csv(DGE, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS4c.csv"),row.names = T)


##################################################### readthru



#### using genes with only 1 PAS

# The log2(RPM of FT sample/RPM of 4sU sample) value is called Stability Score (SS). Two biological replicates were averaged.

library(reshape2)
library(tidyverse)
library(ggrepel)
List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")

plot <- list()
for (cell in c("HepG2_1_7h")) {
  
  PAS_aaa0 <- PAS_aaa <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_APA_strPAS/result/HepG2/pas_17.4.csv"), header = T)
  PAS_aaa0$sumCount <- rowSums(PAS_aaa0[,grepl("_usage$", names(PAS_aaa0))])/length(names(PAS_aaa0)[grepl("_usage$", names(PAS_aaa0))])   ### 20220602 LW, use usage average to determine the most expressed PASs
  
  filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/JTEonly_DGE_mtNormalization/QSF/_DGE_mtNormalization/", pattern="*.csv", full.names=TRUE)[c(10)]
  filenames
  DGE<-List.RS[[4]][c(1,13)]
  names(DGE)[c(2)] <- "RS"
  
  
  
  
  
  PAS_N1 <- subset(PAS_aaa0, region=="3UTR" & sumCount >= 95) %>% count(gene_symbol)
  PAS_aaa1 <- subset(subset(PAS_aaa0, region=="3UTR" & sumCount >= 95), gene_symbol %in% subset(PAS_N1, n=1)$gene_symbol)
  
  
  PAS_N2 <- subset(PAS_aaa0, region=="3UTR" & sumCount >= 10) %>% count(gene_symbol)
  PAS_aaa2 <- subset(subset(PAS_aaa0, region=="3UTR" & sumCount >= 10), gene_symbol %in% subset(PAS_N2, n>=2)$gene_symbol)
  
  
  library(RColorBrewer)
  
  plot[[cell]] <- ggplot() +
    stat_ecdf(data=DGE, aes(x=RS, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol), aes(x=RS, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol), aes(x=RS, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    scale_x_continuous(limits = c(-1, 3)) +
    labs(
      title = paste0("DGE, ",cell,", PAS number is based on 3'UTR",
                     "\nAllvSPA, ", formatC(ks.test(DGE$RS, subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$RS)$p.value, format = "e", digits = 2),
                     "\nAllvAPA, ", formatC(ks.test(DGE$RS, subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$RS)$p.value, format = "e", digits = 2),
                     "\nAPAvSPA, ", formatC(ks.test(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$RS, subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$RS)$p.value, format = "e", digits = 2)
      ),
      x = "delta RS", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("all genes, ",nrow(DGE),", ",round(median(DGE$RS),2)), 
                                  paste0("SPA genes (>=95% usage), ",nrow(subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)),", ",round(median(subset(DGE, gene_symbol %in% PAS_aaa1$gene_symbol)$RS),2)), 
                                  paste0(">=2PAS genes (>=10% usage), ",nrow(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)),", ",round(median(subset(DGE, gene_symbol %in% PAS_aaa2$gene_symbol)$RS),2))),
                       values=c("grey",brewer.pal(n = 8,name = 'Dark2')))
  
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 1)
