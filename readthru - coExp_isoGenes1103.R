
# HeLa, 1 ?M, #8FAADC
# HeLa, 10 ?M, #2F5597
# HepG2, 1 ?M, #F4B183
# HepG2, 10 ?M, #C55A11
# U937, #7030A0
# U937MP, #548235

colorHeLa_1="#8FAADC"
colorHeLa_10="#2F5597"
colorHepG2_1="#F4B183"
colorHepG2_10="#C55A11"
colorU937="#7030A0"
colorU937MP="#548235"





####### get genes whose upstream region 4K does not overlap with upstream genes.

library(reshape)
library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
`%notin%` <- Negate(`%in%`)


isolatedGenes <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/geneSet_isolatedGenes.tbl", header = T, sep = "\t")








########### get commonly expressed (RPM >= 10 in all samples) genes from QSF data

library(tidyverse)
library(RColorBrewer)

List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF.RData")

List.DGE_QSF2 <- lapply(List.DGE_QSF, function(x) { x <- x[c(1,4,5)]; x })
DGE_QSF <- Reduce(function(x, y) merge(x, y), List.DGE_QSF2)
# DGE_QSF2 <- subset(DGE_QSF, DGE_QSF[[2]]>=100 & DGE_QSF[[3]]>=100 & DGE_QSF[[4]]>=100 & DGE_QSF[[5]]>=100 & DGE_QSF[[6]]>=100 & DGE_QSF[[7]]>=100 & DGE_QSF[[8]]>=100 & DGE_QSF[[9]]>=100)
DGE_QSF2 <- subset(DGE_QSF, DGE_QSF[[2]]>=10 & DGE_QSF[[4]]>=10 & DGE_QSF[[6]]>=10 & DGE_QSF[[8]]>=10)

DMSO.common.expressed.genes <- unique(DGE_QSF2$gene_symbol)

# List.DGE_QSF3 <- lapply(List.DGE_QSF, function(x) { x <- x[c(1,10)]; x })
# DGE_QSF3 <- Reduce(function(x, y) merge(x, y), List.DGE_QSF3)
# common.UP.genes <- unique(subset(DGE_QSF3, DGE_QSF3[[2]]=="UP" & DGE_QSF3[[3]]=="UP" & DGE_QSF3[[4]]=="UP" & DGE_QSF3[[5]]=="UP")$gene_symbol)
# common.DN.genes <- unique(subset(DGE_QSF3, DGE_QSF3[[2]]=="DN" & DGE_QSF3[[3]]=="DN" & DGE_QSF3[[4]]=="DN" & DGE_QSF3[[5]]=="DN")$gene_symbol)





List.DGE_QSF4 <- lapply(List.DGE_QSF, function(x) { x$DMSO_common_express <- ifelse(x$gene_symbol %in% DMSO.common.expressed.genes, "yes", "no"); x })
# List.DGE_QSF4 <- lapply(List.DGE_QSF4, function(x) { x$common.UP <- ifelse(x$gene_symbol %in% common.UP.genes, "yes", "no"); x })
# List.DGE_QSF4 <- lapply(List.DGE_QSF4, function(x) { x$common.DN <- ifelse(x$gene_symbol %in% common.DN.genes, "yes", "no"); x })











library(tidyverse)
library(RColorBrewer)

List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")
# List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_4K_DR_4K.RData")

# List.RS2 <- lapply(List.RS, function(x) { x$UR4K_overlap <- ifelse(x$gene_symbol %in% UR4K.olp.genes, "yes", "no"); x })
List.RS2 <- lapply(List.RS, function(x) { x$isolatedGenes <- ifelse(x$gene_symbol %in% isolatedGenes$gene_symbol, "yes", "no"); x })






cells <- c("UP", "U937", "HeLa", "HepG2")

ReportList <- list()

for (cell in cells) {
  
  # cell="UP"
  
  
  
  DGE <- List.DGE_QSF4[[cell]][c(1,4:6,9:11)]
  
  RS <- List.RS2[[cell]][c(1,7:10,13,14)]
  
  
df <- merge(RS, DGE, by="gene_symbol", all=T)
df[is.na(df)] <- "NA"

ReportList[[cell]] <- df


  
}

# library(openxlsx)
# wb = createWorkbook()
# 
# for (n in 1:length(ReportList)) {
#   # n=1
#   name <- names(ReportList)[n]
#   
#   # Write to Excel
#   addWorksheet(wb, name)
#   writeData(wb, sheet = n, ReportList[[n]], rowNames = F, colNames = T)
#   
# }
# 
# saveWorkbook(wb, "./_RS_table/RS&DGE_report.xlsx", overwrite = T)










#### QSF_DGE vs readthru


library(tidyverse)
library(RColorBrewer)

List.RS3 <- lapply(List.RS2, function(x) { x$DMSO_common_express <- ifelse(x$gene_symbol %in% DMSO.common.expressed.genes, "yes", "no"); x })

cor.List.RS <- lapply(List.RS3, function(x) { x <- subset(x, isolatedGenes == "yes" & DMSO_common_express == "yes"); x })
common.RS.genes <- Reduce(function(x, y) intersect(x, y), list(cor.List.RS[[c(1)]][[c(1)]],
                                                               cor.List.RS[[c(2)]][[c(1)]],
                                                               cor.List.RS[[c(3)]][[c(1)]],
                                                               cor.List.RS[[c(4)]][[c(1)]]))
cor.List.RS2 <- lapply(cor.List.RS, function(x) { x <- subset(x, gene_symbol %in% common.RS.genes); x })


###
cor.List.RS3 <- lapply(cor.List.RS2, function(x) { x <- x[c(1,13)]; x })
cor.List.RS4 <- list()
for (sample in names(cor.List.RS3)) {
  
  # sample="HeLa_JTE_01"
  
  df <- cor.List.RS3[[sample]]
  names(df)[c(2)] <- paste0(sample,"_",names(df)[c(2)])
  cor.List.RS4[[sample]] <- df
  
}

DF.List.RS4 <- Reduce(function(x, y) merge(x, y, by="gene_symbol"), cor.List.RS4)
DF.List.RS4$median <- apply(DF.List.RS4[c(2:5)], 1, median)

# write.table(DF.List.RS4, paste0("./RS_common_clean.tbl"), row.names = F, col.names = T, sep = "\t")










cor.DGE_QSF <- lapply(List.DGE_QSF4, function(x) { x <- subset(x, DMSO_common_express == "yes"); x })

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
RSlist <- list()
for (cell in cells) {

  # cell="UP"

  ldf <- cor.DGE_QSF[[cell]]
  QSF_UP <- subset(ldf,ldf[[10]]=="UP")
  QSF_NC <- subset(ldf,ldf[[10]]=="NO")
  QSF_DN <- subset(ldf,ldf[[10]]=="DN")


  RS <- cor.List.RS2[[cell]]

  CDFplot[[cell]] <- ggplot() +
    stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
    stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
    stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(cell, ", QSF.DGE vs delta.RS, GB/DR4K",
                        "\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                        "\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
         x = "Delta RS", y = "Cumulative fraction", color = "") +
    scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)),
                                  paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
                                  paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))

    ),
    values=c("grey","red","blue")) + 
    coord_cartesian(xlim = c(-1,3))
    

}

library(cowplot)
plot_grid(plotlist=CDFplot,ncol = 2)




ggplot() +
  stat_ecdf(data=cor.List.RS2[[c(1)]], aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1.2) +
  stat_ecdf(data=cor.List.RS2[[c(2)]], aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1.2) +
  stat_ecdf(data=cor.List.RS2[[c(3)]], aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1.2) +
  stat_ecdf(data=cor.List.RS2[[c(4)]], aes(x=RS.delta.JTEvDMSO, color = "d"),geom = "step", size=1.2) +
  scale_color_manual(values=c(colorU937MP,colorU937,colorHeLa_10,colorHepG2_10),
    labels = c(paste0(names(cor.List.RS2)[c(1)],", ",nrow(cor.List.RS2[[c(1)]]),", ",round(median(cor.List.RS2[[c(1)]][[c(13)]]),2)),
               paste0(names(cor.List.RS2)[c(2)],", ",nrow(cor.List.RS2[[c(2)]]),", ",round(median(cor.List.RS2[[c(2)]][[c(13)]]),2)),
               paste0(names(cor.List.RS2)[c(3)],", ",nrow(cor.List.RS2[[c(3)]]),", ",round(median(cor.List.RS2[[c(3)]][[c(13)]]),2)),
               paste0(names(cor.List.RS2)[c(4)],", ",nrow(cor.List.RS2[[c(4)]]),", ",round(median(cor.List.RS2[[c(4)]][[c(13)]]),2)))) +
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
  coord_cartesian(xlim = c(-0.5,3.5)) +
  labs(title = paste("delta RS, DR4K/GB, isolatedGenes\nDMSO commonly expressed genes (RPM >=10 in all QSF DMSO samples)",sep = ""),
       x = "delta RS", y = "Cumulative fraction", color = "")

ks.test(cor.List.RS2[[c(1)]]$RS.delta.JTEvDMSO, cor.List.RS2[[c(2)]]$RS.delta.JTEvDMSO)


ggplot() +
  stat_ecdf(data=cor.List.RS2[[c(1)]], aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1.2) +
  stat_ecdf(data=cor.List.RS2[[c(2)]], aes(x=RS.delta.JTEvDMSO, color = "d"),geom = "step", size=1.2) +
  scale_color_manual(values=c(colorU937MP,colorU937),
                     labels = c(paste0(names(cor.List.RS2)[c(1)],", ",nrow(cor.List.RS2[[c(1)]]),", ",round(median(cor.List.RS2[[c(1)]][[c(13)]]),2)),
                                paste0(names(cor.List.RS2)[c(2)],", ",nrow(cor.List.RS2[[c(2)]]),", ",round(median(cor.List.RS2[[c(2)]][[c(13)]]),2)))) +
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
  coord_cartesian(xlim = c(-0.5,3.5)) +
  labs(title = paste0("delta RS, DR4K/GB, isolatedGenes\nDMSO commonly expressed genes (RPM >=10 in all QSF DMSO samples)",
                      "\nHepG2vHeLa, ", formatC(ks.test(cor.List.RS2[[c(1)]]$RS.delta.JTEvDMSO, cor.List.RS2[[c(2)]]$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
       x = "delta RS", y = "Cumulative fraction", color = "")


ggplot() +
  stat_ecdf(data=cor.List.RS2[[c(3)]], aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1.2) +
  stat_ecdf(data=cor.List.RS2[[c(4)]], aes(x=RS.delta.JTEvDMSO, color = "d"),geom = "step", size=1.2) +
  scale_color_manual(values=c(colorHeLa_10,colorHepG2_10),
                     labels = c(paste0(names(cor.List.RS2)[c(3)],", ",nrow(cor.List.RS2[[c(3)]]),", ",round(median(cor.List.RS2[[c(3)]][[c(13)]]),2)),
                                paste0(names(cor.List.RS2)[c(4)],", ",nrow(cor.List.RS2[[c(4)]]),", ",round(median(cor.List.RS2[[c(4)]][[c(13)]]),2)))) +
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
  coord_cartesian(xlim = c(-0.5,3.5)) +
  labs(title = paste0("delta RS, DR4K/GB, isolatedGenes\nDMSO commonly expressed genes (RPM >=10 in all QSF DMSO samples)",
                      "\nHepG2vHeLa, ", formatC(ks.test(cor.List.RS2[[c(3)]]$RS.delta.JTEvDMSO, cor.List.RS2[[c(4)]]$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
       x = "delta RS", y = "Cumulative fraction", color = "")


# write.csv(cor.List.RS2[[c(3)]][c(1,13)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig2e_1.csv"),row.names = F)

# write.csv(cor.List.RS2[[c(4)]][c(1,13)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig2e_2.csv"),row.names = F)


# write.csv(cor.List.RS2[[c(1)]][c(1,13)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig7h_1.csv"),row.names = F)

# write.csv(cor.List.RS2[[c(2)]][c(1,13)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig7h_2.csv"),row.names = F)


