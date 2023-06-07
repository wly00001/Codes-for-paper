workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

#### QSF_DGE vs readthru


library(tidyverse)
library(RColorBrewer)


setList_REDu <- readRDS("./setList_REDu_HeLa_HepG2.RData")
setList_REDi <- readRDS("./setList_REDi_HeLa_HepG2.RData")




List.DGE_QSF <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData"))

List.RS <- readRDS(paste0(workplace,"/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData"))



##################### HepG2

DGE <- List.DGE_QSF[[4]][c(1,6)]
names(DGE) <- c("gene_symbol", "L2FC")
RS <- List.RS[[4]]


LE_HepG2 <- subset(DGE, gene_symbol %in% setList_REDu[[4]]$gene_symbol | gene_symbol %in% setList_REDu[[2]]$gene_symbol)
SH_HepG2 <- subset(DGE, gene_symbol %in% setList_REDu[[7]]$gene_symbol | gene_symbol %in% setList_REDu[[5]]$gene_symbol)

LE_HepG2$APA <- rep("LE",nrow(LE_HepG2))
SH_HepG2$APA <- rep("SH",nrow(SH_HepG2))

APA_HepG2 <- rbind(LE_HepG2,SH_HepG2)
# write.csv(APA_HepG2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5e_1.csv"),row.names = T)



ggplot(APA_HepG2, aes(x=APA, y=L2FC, color=APA, fill=APA)) + 
  geom_boxplot(outlier.shape = NA,coef=0, size=2) + 
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(ylim = c(-0.5,0.5)) +
  scale_color_manual(values=c("red","blue")) +
  scale_fill_manual(values=c("pink1","SteelBlue1")) +
  labs(title = paste0("HepG2_10, QSF.DGE vs REDu",
                      "\nLEvSH, wilcox.test, ", formatC(wilcox.test(LE_HepG2$L2FC, SH_HepG2$L2FC)$p.value, format = "e", digits = 2)),
       x = "3'UTR APA", y = "DGE, log2Ratio", color = "")




LE_HepG2 <- subset(RS, gene_symbol %in% setList_REDu[[4]]$gene_symbol | gene_symbol %in% setList_REDu[[2]]$gene_symbol)
SH_HepG2 <- subset(RS, gene_symbol %in% setList_REDu[[7]]$gene_symbol | gene_symbol %in% setList_REDu[[5]]$gene_symbol)

LE_HepG2$APA <- rep("LE",nrow(LE_HepG2))
SH_HepG2$APA <- rep("SH",nrow(SH_HepG2))

APA_HepG2 <- rbind(LE_HepG2,SH_HepG2)
# write.csv(APA_HepG2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5e_2.csv"),row.names = T)



ggplot(APA_HepG2, aes(x=APA, y=RS.delta.JTEvDMSO, color=APA, fill=APA)) + 
  geom_boxplot(outlier.shape = NA,coef=0, size=2) + 
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(ylim = c(0.75,2.25)) +
  scale_color_manual(values=c("red","blue")) +
  scale_fill_manual(values=c("pink1","SteelBlue1")) +
  labs(title = paste0("HepG2_10, delta RTS vs REDu",
                      "\nLEvSH, wilcox.test, ", formatC(wilcox.test(LE_HepG2$RS.delta.JTEvDMSO, SH_HepG2$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
       x = "3'UTR APA", y = "delta RTS", color = "")









SU_HepG2 <- subset(DGE, gene_symbol %in% setList_REDi[[4]]$gene_symbol | gene_symbol %in% setList_REDi[[2]]$gene_symbol)
AC_HepG2 <- subset(DGE, gene_symbol %in% setList_REDi[[7]]$gene_symbol | gene_symbol %in% setList_REDi[[5]]$gene_symbol)

SU_HepG2$APA <- rep("LE",nrow(SU_HepG2))
AC_HepG2$APA <- rep("SH",nrow(AC_HepG2))

APA_HepG2 <- rbind(SU_HepG2,AC_HepG2)
# write.csv(APA_HepG2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5i_1.csv"),row.names = T)



ggplot(APA_HepG2, aes(x=APA, y=L2FC, color=APA, fill=APA)) + 
  geom_boxplot(outlier.shape = NA,coef=0, size=2) + 
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(ylim = c(-1,0.5)) +
  scale_color_manual(values=c("red","blue")) +
  scale_fill_manual(values=c("pink1","SteelBlue1")) +
  labs(title = paste0("HepG2_10, QSF.DGE vs REDi",
                      "\nLEvSH, wilcox.test, ", formatC(wilcox.test(SU_HepG2$L2FC, AC_HepG2$L2FC)$p.value, format = "e", digits = 2)),
       x = "IPA", y = "DGE, log2Ratio", color = "")



SU_HepG2 <- subset(RS, gene_symbol %in% setList_REDi[[4]]$gene_symbol | gene_symbol %in% setList_REDi[[2]]$gene_symbol)
AC_HepG2 <- subset(RS, gene_symbol %in% setList_REDi[[7]]$gene_symbol | gene_symbol %in% setList_REDi[[5]]$gene_symbol)

SU_HepG2$APA <- rep("LE",nrow(SU_HepG2))
AC_HepG2$APA <- rep("SH",nrow(AC_HepG2))

APA_HepG2 <- rbind(SU_HepG2,AC_HepG2)
# write.csv(APA_HepG2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5i_2.csv"),row.names = T)



ggplot(APA_HepG2, aes(x=APA, y=RS.delta.JTEvDMSO, color=APA, fill=APA)) + 
  geom_boxplot(outlier.shape = NA,coef=0, size=2) + 
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  coord_cartesian(ylim = c(0.75,2.25)) +
  scale_color_manual(values=c("red","blue")) +
  scale_fill_manual(values=c("pink1","SteelBlue1")) +
  labs(title = paste0("HepG2_10, delta RTS vs REDi",
                      "\nLEvSH, wilcox.test, ", formatC(wilcox.test(SU_HepG2$RS.delta.JTEvDMSO, AC_HepG2$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
       x = "IPA", y = "delta RTS", color = "")





# ##################### HeLa
# 
# DGE <- List.DGE_QSF[[2]][c(1,6)]
# names(DGE) <- c("gene_symbol", "L2FC")
# RS <- List.RS[[3]]
# 
# 
# LE_HeLa <- subset(DGE, gene_symbol %in% setList_REDu[[3]]$gene_symbol | gene_symbol %in% setList_REDu[[2]]$gene_symbol)
# SH_HeLa <- subset(DGE, gene_symbol %in% setList_REDu[[6]]$gene_symbol | gene_symbol %in% setList_REDu[[5]]$gene_symbol)
# 
# ggplot() +
#   stat_ecdf(data=DGE, aes(x=L2FC, color = "a"),geom = "step", size=1) +
#   stat_ecdf(data=LE_HeLa, aes(x=L2FC, color = "d"),geom = "step", size=1) +
#   stat_ecdf(data=SH_HeLa, aes(x=L2FC, color = "g"),geom = "step", size=1) +
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         plot.title = element_text(size = 10),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   labs(title = paste0("HeLa_10, QSF.DGE vs REDu",
#                       "\nALLvLE_HeLa, ", formatC(ks.test(DGE$L2FC, LE_HeLa$L2FC)$p.value, format = "e", digits = 2),
#                       "\nALLvSH_HeLa, ", formatC(ks.test(DGE$L2FC, SH_HeLa$L2FC)$p.value, format = "e", digits = 2),
#                       "\nLE_HeLavSH_HeLa, ", formatC(ks.test(LE_HeLa$L2FC, SH_HeLa$L2FC)$p.value, format = "e", digits = 2)),
#        x = "DGE, log2Ratio", y = "Cumulative fraction", color = "") +
#   scale_color_manual(labels = c(paste0("ALL, ",nrow(DGE),", ",round(median(DGE$L2FC),2)), 
#                                 paste0("LE_HeLa, ",nrow(LE_HeLa),", ",round(median(LE_HeLa$L2FC),2)), 
#                                 paste0("SH_HeLa, ",nrow(SH_HeLa),", ",round(median(SH_HeLa$L2FC),2))), 
#                      values=c("grey",brewer.pal(n = 12,name = 'Paired')[c(5,1)])) + scale_x_continuous(limits = c(-2, 1))
# 
# 
# LE_HeLa <- subset(RS, gene_symbol %in% setList_REDu[[3]]$gene_symbol | gene_symbol %in% setList_REDu[[2]]$gene_symbol)
# SH_HeLa <- subset(RS, gene_symbol %in% setList_REDu[[6]]$gene_symbol | gene_symbol %in% setList_REDu[[5]]$gene_symbol)
# 
# ggplot() +
#   stat_ecdf(data=RS, aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
#   stat_ecdf(data=LE_HeLa, aes(x=RS.delta.JTEvDMSO, color = "d"),geom = "step", size=1) +
#   stat_ecdf(data=SH_HeLa, aes(x=RS.delta.JTEvDMSO, color = "g"),geom = "step", size=1) +
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         plot.title = element_text(size = 10),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   labs(title = paste0("HeLa_1, RS vs REDu",
#                       "\nALLvLE_HeLa, ", formatC(ks.test(RS$RS.delta.JTEvDMSO, LE_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
#                       "\nALLvSH_HeLa, ", formatC(ks.test(RS$RS.delta.JTEvDMSO, SH_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
#                       "\nLE_HeLavSH_HeLa, ", formatC(ks.test(LE_HeLa$RS.delta.JTEvDMSO, SH_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
#        x = "Readthrough score", y = "Cumulative fraction", color = "") + 
#   scale_x_continuous(limits = c(-1, 3)) +
#   scale_color_manual(labels = c(paste0("ALL, ",nrow(RS),", ",round(median(RS$RS.delta.JTEvDMSO),2)), 
#                                 paste0("LE_HeLa, ",nrow(LE_HeLa),", ",round(median(LE_HeLa$RS.delta.JTEvDMSO),2)), 
#                                 paste0("SH_HeLa, ",nrow(SH_HeLa),", ",round(median(SH_HeLa$RS.delta.JTEvDMSO),2))), 
#                      values=c("grey",brewer.pal(n = 12,name = 'Paired')[c(5,1)]))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# SU_HeLa <- subset(DGE, gene_symbol %in% setList_REDi[[3]]$gene_symbol | gene_symbol %in% setList_REDi[[2]]$gene_symbol)
# AC_HeLa <- subset(DGE, gene_symbol %in% setList_REDi[[6]]$gene_symbol | gene_symbol %in% setList_REDi[[5]]$gene_symbol)
# 
# ggplot() +
#   stat_ecdf(data=DGE, aes(x=L2FC, color = "a"),geom = "step", size=1) +
#   stat_ecdf(data=SU_HeLa, aes(x=L2FC, color = "d"),geom = "step", size=1) +
#   stat_ecdf(data=AC_HeLa, aes(x=L2FC, color = "g"),geom = "step", size=1) +
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         plot.title = element_text(size = 10),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   labs(title = paste0("HeLa_10, QSF.DGE vs REDi",
#                       "\nALLvSU_HeLa, ", formatC(ks.test(DGE$L2FC, SU_HeLa$L2FC)$p.value, format = "e", digits = 2),
#                       "\nALLvAC_HeLa, ", formatC(ks.test(DGE$L2FC, AC_HeLa$L2FC)$p.value, format = "e", digits = 2),
#                       "\nSU_HeLavAC_HeLa, ", formatC(ks.test(SU_HeLa$L2FC, AC_HeLa$L2FC)$p.value, format = "e", digits = 2)),
#        x = "DGE, log2Ratio", y = "Cumulative fraction", color = "") +
#   scale_color_manual(labels = c(paste0("ALL, ",nrow(DGE),", ",round(median(DGE$L2FC),2)), 
#                                 paste0("SU_HeLa, ",nrow(SU_HeLa),", ",round(median(SU_HeLa$L2FC),2)), 
#                                 paste0("AC_HeLa, ",nrow(AC_HeLa),", ",round(median(AC_HeLa$L2FC),2))), 
#                      values=c("grey",brewer.pal(n = 12,name = 'Paired')[c(5,1)])) + scale_x_continuous(limits = c(-2, 1))
# 
# 
# 
# SU_HeLa <- subset(RS, gene_symbol %in% setList_REDi[[3]]$gene_symbol | gene_symbol %in% setList_REDi[[2]]$gene_symbol)
# AC_HeLa <- subset(RS, gene_symbol %in% setList_REDi[[6]]$gene_symbol | gene_symbol %in% setList_REDi[[5]]$gene_symbol)
# 
# ggplot() +
#   stat_ecdf(data=RS, aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
#   stat_ecdf(data=SU_HeLa, aes(x=RS.delta.JTEvDMSO, color = "d"),geom = "step", size=1) +
#   stat_ecdf(data=AC_HeLa, aes(x=RS.delta.JTEvDMSO, color = "g"),geom = "step", size=1) +
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         plot.title = element_text(size = 10),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   labs(title = paste0("HeLa_1, RS vs REDi",
#                       "\nALLvSU_HeLa, ", formatC(ks.test(RS$RS.delta.JTEvDMSO, SU_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
#                       "\nALLvAC_HeLa, ", formatC(ks.test(RS$RS.delta.JTEvDMSO, AC_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
#                       "\nSU_HeLavAC_HeLa, ", formatC(ks.test(SU_HeLa$RS.delta.JTEvDMSO, AC_HeLa$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
#        x = "Readthrough score", y = "Cumulative fraction", color = "") + 
#   scale_x_continuous(limits = c(-1, 3)) +
#   scale_color_manual(labels = c(paste0("ALL, ",nrow(RS),", ",round(median(RS$RS.delta.JTEvDMSO),2)), 
#                                 paste0("SU_HeLa, ",nrow(SU_HeLa),", ",round(median(SU_HeLa$RS.delta.JTEvDMSO),2)), 
#                                 paste0("AC_HeLa, ",nrow(AC_HeLa),", ",round(median(AC_HeLa$RS.delta.JTEvDMSO),2))), 
#                      values=c("grey",brewer.pal(n = 12,name = 'Paired')[c(5,1)]))
