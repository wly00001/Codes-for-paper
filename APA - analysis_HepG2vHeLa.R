
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


########### HepG2 v HeLa, U937MP v U937

library(reshape2)
library(tidyverse)
library(ggrepel)

List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDu.RData"))[c(5)]
List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDi.RData"))[c(5)]

# List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_U937_U937P/List.REDu.RData"))[c(3)]
# List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_U937_U937P/List.REDi.RData"))[c(3)]



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
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
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
plot_grid(plotlist=plot_REDi,ncol = 1)


# write.csv(List.REDu[[1]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig6b.csv"),row.names = T)

# write.csv(List.REDi[[1]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS5a.csv"),row.names = T)








########### HepG2 v HeLa, U937MP v U937, RE CDF plots

library(reshape2)
library(tidyverse)
library(ggrepel)

List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDu.RData"))[c(5)]
List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDi.RData"))[c(5)]

# List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_U937_U937P/List.REDu.RData"))[c(3)]
# List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_U937_U937P/List.REDi.RData"))[c(3)]



List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x<-na.omit(x); x })
List.REDi <- lapply(List.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDi <- lapply(List.REDi, function(x) { x<-na.omit(x); x })






ggplot() +
  stat_ecdf(data=List.REDu[[1]], aes(x=RE_test, color = "a"),geom = "step", size=1.2) +
  stat_ecdf(data=List.REDu[[1]], aes(x=RE_ctrl, color = "b"),geom = "step", size=1.2) +
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
  labs(title = paste0("APA, RE=dPAS/pPAS",
                      "\nHepG2vHeLa, ", formatC(ks.test(List.REDu[[1]]$RE_test, List.REDu[[1]]$RE_ctrl)$p.value, format = "e", digits = 2)),
       x = "RE, log2(dPAS/pPAS)", y = "Cumulative fraction", color = "") +
  scale_color_manual(labels = c(paste0("HepG2 DMSO, ",nrow(List.REDu[[1]]),", ",round(median(List.REDu[[1]]$RE_test),2)),
                                paste0("HeLa DMSO, ",nrow(List.REDu[[1]]),", ",round(median(List.REDu[[1]]$RE_ctrl),2))), 
                     values=c(colorHepG2_10,colorHeLa_10)) + 
  coord_cartesian(xlim = c(-4, 3))











library(reshape2)
library(tidyverse)
library(ggrepel)

List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDi.RData"))[c(4)]

PAS_aaa1 <- PAS_aaa <- read.csv(paste0("./APA_HeLa_HepG2/pas_17.4.csv"), header = T)

MRPL12 <- subset(PAS_aaa1, gene_symbol == "MRPL12")



contingency_table <- matrix(c(127,sum(PAS_aaa1$HepG2_DMSO_count)-127,
                              149,sum(PAS_aaa1$HepG2_JTE10_count)-149),nrow=2)
fisher.test(contingency_table)$p.value


contingency_table <- matrix(c(586,sum(PAS_aaa1$HepG2_DMSO_count)-586,
                              253,sum(PAS_aaa1$HepG2_JTE10_count)-253),nrow=2)
fisher.test(contingency_table)$p.value

