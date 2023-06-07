
workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"




####### APA CDF

library(tidyverse)
library(RColorBrewer)


# List.REDu <- readRDS(paste0("./result/FIP1_JTE_PASS/List.REDu.RData"))
# List.REDi <- readRDS(paste0("./result/FIP1_JTE_PASS/List.REDi.RData"))
List.REDu <- readRDS(paste0(workplace,"/linux/lwang/21029-09/R2/R2_uniqueMapping/result/FIP1_JTE_0327/List.REDu.RData"))
List.REDi <- readRDS(paste0(workplace,"/linux/lwang/21029-09/R2/R2_uniqueMapping/result/FIP1_JTE_0327/List.REDi_1vAll.RData"))

names(List.REDu) <- gsub("YC_103_pA_1_1","24DnJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_2","24DnJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_3","24DpJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_4","24DpJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_5","48DnJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_6","48DnJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_7","48DpJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_8","48DpJp",names(List.REDu))

names(List.REDi) <- gsub("YC_103_pA_1_1","24DnJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_2","24DnJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_3","24DpJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_4","24DpJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_5","48DnJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_6","48DnJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_7","48DpJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_8","48DpJp",names(List.REDi))


List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x<-na.omit(x); x })
List.REDi <- lapply(List.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDi <- lapply(List.REDi, function(x) { x<-na.omit(x); x })




List.REDu2 <- list()
List.REDi2 <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDu2[[sample]] <- df
  
}

for (sample in names(List.REDi)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDi[[sample]]
  names(df)[c(4:13)] <- paste0(names(df)[c(4:13)],".",sample)
  List.REDi2[[sample]] <- df
  
}






################################################################################
################################################################################
################################################################################

# APA plot


ratio_tbl_REDu <- data.frame(comparison=character(),
                             FoSpm_TPA_LE=numeric(),
                             FoSpm_TPA_SH=numeric(),
                             TPA_Ratio=numeric(),
                             stringsAsFactors=FALSE)

plot_REDu <- list()
for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  REDu_aaa <-List.REDu[[sample]]
  plot_REDu[[sample]] <- ggplot() +
    geom_point(data=subset(REDu_aaa, regu.adj == "NO"), aes(x=as.numeric(RE_pPAS), y=as.numeric(RE_dPAS), color="NO"), size = 0.1, alpha=0.3)+
    geom_point(data=subset(REDu_aaa, regu.adj == "lengthened"), aes(x=as.numeric(RE_pPAS), y=as.numeric(RE_dPAS), color="lengthened"), size = 1.5, alpha=0.3)+
    geom_point(data=subset(REDu_aaa, regu.adj == "shortened"), aes(x=as.numeric(RE_pPAS), y=as.numeric(RE_dPAS), color="shortened"), size = 1.5, alpha=0.3)+
    scale_color_manual(values=c("red", "grey", "blue")) +
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          legend.position = "none",
          plot.title = element_text(size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    scale_y_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    # annotate('text', x=-2, y=2, label=nrow(subset(REDu_aaa, regu.adj == "lengthened")),color = "red") +
    # annotate('text', x=2, y=-2, label=nrow(subset(REDu_aaa, regu.adj == "shortened")),color = "blue") +
    # geom_label_repel(data=subset(REDu_aaa, gene %in% c("TIMP2")),
    #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
    labs(title = paste("REDu, ",sample,
                       "\nred=",table(REDu_aaa$regu.adj)[1],
                       ", gray=",table(REDu_aaa$regu.adj)[2],
                       ", blue=",table(REDu_aaa$regu.adj)[3],
                       ", ratio= ",ifelse(table(REDu_aaa$regu.adj)[1]>table(REDu_aaa$regu.adj)[3],round(table(REDu_aaa$regu.adj)[1]/table(REDu_aaa$regu.adj)[3],1),(-1)*round(table(REDu_aaa$regu.adj)[3]/table(REDu_aaa$regu.adj)[1],1)),
                       sep = ""),
         x = "Log2Ratio of pPAS isoforms", y = "Log2Ratio of dPAS isoforms", color = "")
  
  ratio_tbl_REDu[nrow(ratio_tbl_REDu) + 1,] = c(sample,table(REDu_aaa$regu.adj)[1],table(REDu_aaa$regu.adj)[3],ifelse(table(REDu_aaa$regu.adj)[1]>table(REDu_aaa$regu.adj)[3],round(table(REDu_aaa$regu.adj)[1]/table(REDu_aaa$regu.adj)[3],1),(-1)*round(table(REDu_aaa$regu.adj)[3]/table(REDu_aaa$regu.adj)[1],1)))
  
}



ratio_tbl_REDi <- data.frame(comparison=character(),
                             FoSpm_IPA_SU=numeric(),
                             FoSpm_IPA_AC=numeric(),
                             IPA_Ratio=numeric(),
                             stringsAsFactors=FALSE)

plot_REDi <- list()
for (sample in names(List.REDi)) {
  
  # sample="HeLa_JTE_01"
  
  REDi_aaa <-List.REDi[[sample]]
  plot_REDi[[sample]] <- ggplot() +
    geom_point(data=subset(REDi_aaa, regu.adj == "NO"), aes(x=as.numeric(RE_IPA), y=as.numeric(RE_TPA), color="NO"), size = 0.1, alpha=0.3)+
    geom_point(data=subset(REDi_aaa, regu.adj == "lengthened"), aes(x=as.numeric(RE_IPA), y=as.numeric(RE_TPA), color="lengthened"), size = 1.5, alpha=0.3)+
    geom_point(data=subset(REDi_aaa, regu.adj == "shortened"), aes(x=as.numeric(RE_IPA), y=as.numeric(RE_TPA), color="shortened"), size = 1.5, alpha=0.3)+
    scale_color_manual(values=c("red", "grey", "blue")) +
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          legend.position = "none",
          plot.title = element_text(size = 8),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    scale_y_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    # annotate('text', x=-2, y=2, label=nrow(subset(REDi_aaa, regu.adj == "lengthened")),color = "red") +
    # annotate('text', x=2, y=-2, label=nrow(subset(REDi_aaa, regu.adj == "shortened")),color = "blue") +
    # geom_label_repel(data=subset(REDi_aaa, gene %in% c("TIMP2")),
    #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
    labs(title = paste("REDi, ",sample,
                       "\nred=",table(REDi_aaa$regu.adj)[1],
                       ", gray=",table(REDi_aaa$regu.adj)[2],
                       ", blue=",table(REDi_aaa$regu.adj)[3],
                       ", ratio= ",ifelse(table(REDi_aaa$regu.adj)[1]>table(REDi_aaa$regu.adj)[3],round(table(REDi_aaa$regu.adj)[1]/table(REDi_aaa$regu.adj)[3],1),(-1)*round(table(REDi_aaa$regu.adj)[3]/table(REDi_aaa$regu.adj)[1],1)),
                       sep = ""),
         x = "Log2Ratio of IPA isoforms", y = "Log2Ratio of TPA isoforms", color = "")
  
  ratio_tbl_REDi[nrow(ratio_tbl_REDi) + 1,] = c(sample,table(REDi_aaa$regu.adj)[1],table(REDi_aaa$regu.adj)[3],ifelse(table(REDi_aaa$regu.adj)[1]>table(REDi_aaa$regu.adj)[3],round(table(REDi_aaa$regu.adj)[1]/table(REDi_aaa$regu.adj)[3],1),(-1)*round(table(REDi_aaa$regu.adj)[3]/table(REDi_aaa$regu.adj)[1],1)))
  
}

library(cowplot)
plot_grid(plotlist=c(plot_REDu[c(7)],plot_REDi[c(7)]),ncol = 1)

# write.csv(List.REDu[[7]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig9g.csv"),row.names = T)
# write.csv(List.REDi[[7]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS9d.csv"),row.names = T)













ldf.REDu <- lapply(List.REDu, function(x) { x <- x[c(1,9,13,14)]; x })
ldf.REDu <- lapply(ldf.REDu, function(x) { names(x) <- c("gene_symbol","RED","pval.fisher.adj","regu.adj"); x })
# ldf.REDi <- lapply(List.REDi, function(x) { x <- x[c(1,9,13,14)]; x })
ldf.REDi <- lapply(List.REDi, function(x) { x <- x[c(1,6,10,11)]; x })
ldf.REDi <- lapply(ldf.REDi, function(x) { names(x) <- c("gene_symbol","RED","pval.fisher.adj","regu.adj"); x })


ggplot() +
  stat_ecdf(data=ldf.REDu[[c(6)]], aes(x=RED, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=ldf.REDu[[c(10)]], aes(x=RED, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(brewer.pal(n = 12,name = 'Paired')),
                     labels = c(paste0(names(ldf.REDu)[c(6)],", ",nrow(ldf.REDu[[c(6)]]),", ",round(median(ldf.REDu[[c(6)]][[c(2)]]),2)),
                                paste0(names(ldf.REDu)[c(10)],", ",nrow(ldf.REDu[[c(10)]]),", ",round(median(ldf.REDu[[c(10)]][[c(2)]]),2)))) +
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
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(title = paste0("REDu, PolyA_DB.v4s",
                      "\nKS test, ", formatC(ks.test(ldf.REDu[[c(6)]]$RED, ldf.REDu[[c(10)]]$RED)$p.value, format = "e", digits = 2)),
       x = "REDu", y = "Cumulative fraction", color = "")




ggplot() +
  stat_ecdf(data=ldf.REDi[[c(6)]], aes(x=RED, color = "c"),geom = "step", size=1) +
  stat_ecdf(data=ldf.REDi[[c(10)]], aes(x=RED, color = "d"),geom = "step", size=1) +
  scale_color_manual(values=c(brewer.pal(n = 12,name = 'Paired')),
                     labels = c(paste0(names(ldf.REDi)[c(6)],", ",nrow(ldf.REDi[[c(6)]]),", ",round(median(ldf.REDi[[c(6)]][[c(2)]]),2)),
                                paste0(names(ldf.REDi)[c(10)],", ",nrow(ldf.REDi[[c(10)]]),", ",round(median(ldf.REDi[[c(10)]][[c(2)]]),2)))) +
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
  coord_cartesian(xlim = c(-3, 2)) +
  labs(title = paste0("REDi, PolyA_DB.v4s",
                      "\nKS test, ", formatC(ks.test(ldf.REDi[[c(6)]]$RED, ldf.REDi[[c(10)]]$RED)$p.value, format = "e", digits = 2)),
       x = "REDi", y = "Cumulative fraction", color = "")









#######  APA colored UpSet plot

library(reshape2)
library(tidyverse)
library(RColorBrewer)

List.REDu <- readRDS(paste0(workplace,"/linux/lwang/21029-09/R2/R2_uniqueMapping/result/FIP1_JTE_0327/List.REDu.RData"))
List.REDi <- readRDS(paste0(workplace,"/linux/lwang/21029-09/R2/R2_uniqueMapping/result/FIP1_JTE_0327/List.REDi_1vAll.RData"))

names(List.REDu) <- gsub("YC_103_pA_1_1","24DnJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_2","24DnJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_3","24DpJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_1_4","24DpJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_5","48DnJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_6","48DnJp",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_7","48DpJn",names(List.REDu))
names(List.REDu) <- gsub("YC_103_pA_2_8","48DpJp",names(List.REDu))

names(List.REDi) <- gsub("YC_103_pA_1_1","24DnJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_2","24DnJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_3","24DpJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_1_4","24DpJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_5","48DnJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_6","48DnJp",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_7","48DpJn",names(List.REDi))
names(List.REDi) <- gsub("YC_103_pA_2_8","48DpJp",names(List.REDi))


List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x<-na.omit(x); x })
List.REDi <- lapply(List.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDi <- lapply(List.REDi, function(x) { x<-na.omit(x); x })




List.REDu2 <- list()
List.REDi2 <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDu2[[sample]] <- df
  
}

for (sample in names(List.REDi)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDi[[sample]]
  names(df)[c(4:13)] <- paste0(names(df)[c(4:13)],".",sample)
  List.REDi2[[sample]] <- df
  
}




List.REDu.8h3 <- lapply(List.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3[c(6,10)])
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO")

library(UpSetR)
list.Dox <- list(sig.LE_Dn = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="lengthened")$gene_symbol,
                          sig.LE_Dp = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="lengthened")$gene_symbol,
                          nonsig.LE_Dn = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]>0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.LE_Dp = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]>0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol,
                          sig.SH_Dn = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="shortened")$gene_symbol,
                          sig.SH_Dp = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="shortened")$gene_symbol,
                          nonsig.SH_Dn = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]<0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.SH_Dp = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]<0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol)
upset(fromList(list.Dox), nsets = 8, nintersects = NA, sets = names(list.Dox), keep.order = T, order.by = "freq", text.scale = 1.5)

# write.csv(heatmapTbl.REDu2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig9j.csv"),row.names = T)




List.REDi.8h3 <- lapply(List.REDi2, function(x) { x <- x[c(1,8,13)]; x })
heatmapTbl.REDi <- Reduce(function(x, y) merge(x, y), List.REDi.8h3[c(6,10)])
heatmapTbl.REDi2 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO" | heatmapTbl.REDi[[5]]!="NO")

library(UpSetR)
list.Dox <- list(sig.SU_Dn = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="lengthened")$gene_symbol,
                          sig.SU_Dp = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="lengthened")$gene_symbol,
                          nonsig.SU_Dn = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]>0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.SU_Dp = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]>0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol,
                          sig.AC_Dn = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="shortened")$gene_symbol,
                          sig.AC_Dp = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="shortened")$gene_symbol,
                          nonsig.AC_Dn = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]<0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.AC_Dp = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]<0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol)
upset(fromList(list.Dox), nsets = 8, nintersects = NA, sets = names(list.Dox), keep.order = T, order.by = "freq", text.scale = 1.5)

# write.csv(heatmapTbl.REDi2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS9e.csv"),row.names = T)

