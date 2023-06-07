
#### QSF_DGE vs readthru


library(tidyverse)
library(RColorBrewer)

List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData")
# List.DGE_ChrRNA <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_ChrRNA.RData")

List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")
# List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_4K_DR_4K.RData")

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
RSlist <- list()
for (cell in cells) {
  
  # cell="HepG2"
  
  RS <- List.RS[[cell]]
  
  if (cell == "HeLa") {
    
    for (c in c("HeLa_1", "HeLa_10")){
      
      ldf <- List.DGE_QSF[[c]]
      QSF_UP <- subset(ldf,ldf[[10]]=="UP")
      QSF_NC <- subset(ldf,ldf[[10]]=="NO")
      QSF_DN <- subset(ldf,ldf[[10]]=="DN")
      
      
      
      CDFplot[[c]] <- ggplot() +
        # stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
        labs(title = paste0(c, ", QSF.DGE vs delta.RS, GB/DR4K",
                            #"\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            #"\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            "\nUPvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
             x = "Delta RS", y = "Cumulative fraction", color = "") +
        scale_color_manual(labels = c(#paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
          paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
          paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
          
        ), 
        # values=c("grey","red","blue")) + 
        values=c("red","blue")) + 
        coord_cartesian(xlim = c(-1, 3))
      
      # write.csv(subset(RS, gene_symbol %in% QSF_UP$gene_symbol|gene_symbol %in% QSF_DN$gene_symbol), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig2f_",c,".csv"),row.names = F)
      
      
      
      
    }
    
    
    
  } else if (cell == "HepG2") {
    
    for (c in c("HepG2_1", "HepG2_10")){
      
      ldf <- List.DGE_QSF[[c]]
      QSF_UP <- subset(ldf,ldf[[10]]=="UP")
      QSF_NC <- subset(ldf,ldf[[10]]=="NO")
      QSF_DN <- subset(ldf,ldf[[10]]=="DN")
      
      
      
      CDFplot[[c]] <- ggplot() +
        # stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
        labs(title = paste0(c, ", QSF.DGE vs delta.RS, GB/DR4K",
                            #"\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            #"\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            "\nUPvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
             x = "Delta RS", y = "Cumulative fraction", color = "") +
        scale_color_manual(labels = c(#paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
          paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
          paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
          
        ), 
        # values=c("grey","red","blue")) + 
        values=c("red","blue")) + 
        coord_cartesian(xlim = c(-1, 3))
      
      
      # write.csv(subset(RS, gene_symbol %in% QSF_UP$gene_symbol|gene_symbol %in% QSF_DN$gene_symbol), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig2f_",c,".csv"),row.names = F)
      
      
      
    }
    
    
    
  } else {
    
    ldf <- List.DGE_QSF[[cell]]
    QSF_UP <- subset(ldf,ldf[[10]]=="UP")
    QSF_NC <- subset(ldf,ldf[[10]]=="NO")
    QSF_DN <- subset(ldf,ldf[[10]]=="DN")
    
    
    CDFplot[[cell]] <- ggplot() +
      # stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
      stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
      stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
      labs(title = paste0(cell, ", QSF.DGE vs delta.RS, GB/DR4K",
                          #"\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                          #"\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                          "\nUPvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
           x = "Delta RS", y = "Cumulative fraction", color = "") +
      scale_color_manual(labels = c(#paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
        paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
        paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
        
      ), 
      # values=c("grey","red","blue")) + 
      values=c("red","blue")) + 
      coord_cartesian(xlim = c(-1, 3))
    
    
    
    
  }
  
  
}

library(cowplot)
# plot_grid(plotlist=CDFplot,ncol = 2)
plot_grid(plotlist=CDFplot[c(4)],ncol = 1) #for paper
plot_grid(plotlist=CDFplot[c(6)],ncol = 1) #for paper

plot_grid(plotlist=CDFplot[c(1)],ncol = 1) #for paper
plot_grid(plotlist=CDFplot[c(2)],ncol = 1) #for paper

















#### QSF_DGE vs readthru for isolated genes


library(tidyverse)
library(RColorBrewer)


isolatedGenes <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/geneSet_isolatedGenes.tbl", header = T, sep = "\t")



List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData")
# List.DGE_ChrRNA <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_ChrRNA.RData")

List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")
# List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_4K_DR_4K.RData")

cells <- c("UP", "U937", "HeLa", "HepG2")

CDFplot <- list()
RSlist <- list()
for (cell in cells) {
  
  # cell="HeLa"
  
  RS <- List.RS[[cell]]
  
  if (cell == "HeLa") {
    
    for (c in c("HeLa_1", "HeLa_10")){
      
      # c="HeLa_1"
      
      ldf <- List.DGE_QSF[[c]]
      ldf <- subset(ldf, gene_symbol %in% isolatedGenes$gene_symbol)
      
      QSF_UP <- subset(ldf,ldf[[10]]=="UP")
      QSF_NC <- subset(ldf,ldf[[10]]=="NO")
      QSF_DN <- subset(ldf,ldf[[10]]=="DN")
      
      
      
      CDFplot[[c]] <- ggplot() +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
        labs(title = paste0(c, ", QSF.DGE vs delta.RS, GB/DR4K",
                            "\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            "\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
             x = "Delta RS", y = "Cumulative fraction", color = "") +
        scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
                                      paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
                                      paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
                                      
        ), 
        values=c("grey","red","blue")) + coord_cartesian(xlim = c(-2, 4))
      
      
      
      
      
    }
    
    
    
  } else if (cell == "HepG2") {
    
    for (c in c("HepG2_1", "HepG2_10")){
      
      ldf <- List.DGE_QSF[[c]]
      ldf <- subset(ldf, gene_symbol %in% isolatedGenes$gene_symbol)
      
      QSF_UP <- subset(ldf,ldf[[10]]=="UP")
      QSF_NC <- subset(ldf,ldf[[10]]=="NO")
      QSF_DN <- subset(ldf,ldf[[10]]=="DN")
      
      
      
      CDFplot[[c]] <- ggplot() +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
        stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
        labs(title = paste0(c, ", QSF.DGE vs delta.RS, GB/DR4K",
                            "\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                            "\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
             x = "Delta RS", y = "Cumulative fraction", color = "") +
        scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
                                      paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
                                      paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
                                      
        ), 
        values=c("grey","red","blue")) + coord_cartesian(xlim = c(-1, 3))
      
      
      
      
      
    }
    
    
    
  } else {
    
    ldf <- List.DGE_QSF[[cell]]
    ldf <- subset(ldf, gene_symbol %in% isolatedGenes$gene_symbol)
    
    QSF_UP <- subset(ldf,ldf[[10]]=="UP")
    QSF_NC <- subset(ldf,ldf[[10]]=="NO")
    QSF_DN <- subset(ldf,ldf[[10]]=="DN")
    
    
    CDFplot[[cell]] <- ggplot() +
      stat_ecdf(data=subset(RS, gene_symbol %in% QSF_NC$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
      stat_ecdf(data=subset(RS, gene_symbol %in% QSF_UP$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
      stat_ecdf(data=subset(RS, gene_symbol %in% QSF_DN$gene_symbol), aes(x=RS.delta.JTEvDMSO, color = "c"),geom = "step", size=1) +
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
      labs(title = paste0(cell, ", QSF.DGE vs delta.RS, GB/DR4K",
                          "\nNCvUP, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2),
                          "\nNCvDN, ", formatC(ks.test(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO, subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO)$p.value, format = "e", digits = 2)),
           x = "Delta RS", y = "Cumulative fraction", color = "") +
      scale_color_manual(labels = c(paste0("NC genes, ",nrow(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_NC$gene_symbol)$RS.delta.JTEvDMSO),2)), 
                                    paste0("UP genes, ",nrow(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_UP$gene_symbol)$RS.delta.JTEvDMSO),2)),
                                    paste0("DN genes, ",nrow(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)),", ",round(median(subset(RS, gene_symbol %in% QSF_DN$gene_symbol)$RS.delta.JTEvDMSO),2))
                                    
      ), 
      values=c("grey","red","blue")) + coord_cartesian(xlim = c(-2, 4))
    
    
    
  }
  
  
}

library(cowplot)
plot_grid(plotlist=CDFplot,ncol = 2)