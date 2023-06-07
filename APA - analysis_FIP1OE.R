workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"


library(reshape2)
library(tidyverse)
library(ggrepel)

List.REDu <- readRDS(paste0(workplace,"/linux/lwang/archive/FIP1L1/FIP1L1_OE/result/OE_0406/List.REDu.RData"))
List.REDi <- readRDS(paste0(workplace,"/linux/lwang/archive/FIP1L1/FIP1L1_OE/result/OE_0406/List.REDi_1vAll.RData"))


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
plot_grid(plotlist=c(plot_REDu,plot_REDi),ncol = length(plot_REDu))

# write.csv(List.REDu[[1]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig8b.csv"),row.names = T)
# write.csv(List.REDi[[1]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig8c.csv"),row.names = T)
