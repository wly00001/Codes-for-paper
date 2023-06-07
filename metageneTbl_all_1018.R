workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

library(reshape)
library(tidyverse)
`%notin%` <- Negate(`%in%`)
CMA <- function(x, n = 1){stats::filter(x, rep(1 / n, n), sides = 2)}

plot <- list()




metaTbl.plus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/_metagene/U937PMA_maxSize.plus.tab"), header = T, sep = "\t")
metaTbl.minus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/_metagene/U937PMA_maxSize.minus.tab"), header = T, sep = "\t")
names(metaTbl.plus) <- gsub(".plus","",names(metaTbl.plus))
names(metaTbl.minus) <- gsub(".minus","",names(metaTbl.minus))

metaTbl <- rbind(metaTbl.plus,metaTbl.minus)

metaTbl.MOCK <- metaTbl[,grep("MOCK",names(metaTbl))]
metaTbl.JTE <- metaTbl[,grep("JTE",names(metaTbl))]


metaTbl.MOCK$rsum <- apply(metaTbl.MOCK, 1, sum,na.rm = T)
metaTbl.MOCK <- metaTbl.MOCK[order(metaTbl.MOCK$rsum,decreasing=T),]
N=round(0.05*nrow(metaTbl.MOCK))
metaTbl.MOCK <- metaTbl.MOCK[-c(1:N,(nrow(metaTbl.MOCK)-N):nrow(metaTbl.MOCK)), ]
# metaTbl.MOCK[c(1:220)] <- metaTbl.MOCK[c(1:220)]/metaTbl.MOCK$rsum
metaTbl.MOCK$rsum <- NULL

metaTbl.JTE$rsum <- apply(metaTbl.JTE, 1, sum,na.rm = T)
metaTbl.JTE <- metaTbl.JTE[order(metaTbl.JTE$rsum,decreasing=T),]
N=round(0.05*nrow(metaTbl.JTE))
metaTbl.JTE <- metaTbl.JTE[-c(1:N,(nrow(metaTbl.JTE)-N):nrow(metaTbl.JTE)), ]
# metaTbl.JTE[c(1:220)] <- metaTbl.JTE[c(1:220)]/metaTbl.JTE$rsum
metaTbl.JTE$rsum <- NULL






plotTbl <- data.frame("DMSO"=apply(metaTbl.MOCK, 2, mean, na.rm = T),"JTE"=apply(metaTbl.JTE, 2, mean, na.rm = T))
plotTbl$pos <- c(1:220)
# plotTbl$DMSO <- as.numeric(CMA(plotTbl$DMSO))
# plotTbl$JTE <- as.numeric(CMA(plotTbl$JTE))
plotTbl <- melt(plotTbl,id.vars = "pos")
names(plotTbl) <- c("pos", "sample", "mean.rpm")
# write.csv(plotTbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/whole_U937MP.csv"),row.names = F)

plot[[paste0("U937MP",".whole")]] <- ggplot(plotTbl, aes(x=pos, y=mean.rpm, group=sample, color=sample)) +
  geom_line(size=1)+
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  scale_x_continuous(limits = c(10,210),breaks=seq(10,210,50),labels=c("-5kb","TSS","","TES","5kb")) +
  scale_color_manual(values=c("grey","red")) + 
  coord_cartesian(ylim = c(0.05, 0.22)) +
  scale_y_continuous(breaks = seq(0.05, 0.2, by=0.025)) +
  labs(title = "U937MP")







metaTbl.plus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/_metagene/TES5K.U937PMA_maxSize.plus.tab"), header = T, sep = "\t")
metaTbl.minus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/_metagene/TES5K.U937PMA_maxSize.minus.tab"), header = T, sep = "\t")
names(metaTbl.plus) <- gsub(".plus","",names(metaTbl.plus))
names(metaTbl.minus) <- gsub(".minus","",names(metaTbl.minus))

metaTbl <- rbind(metaTbl.plus,metaTbl.minus)

metaTbl.MOCK <- metaTbl[,grep("MOCK",names(metaTbl))]
metaTbl.JTE <- metaTbl[,grep("JTE",names(metaTbl))]


metaTbl.MOCK$rsum <- apply(metaTbl.MOCK, 1, sum,na.rm = T)
metaTbl.MOCK <- metaTbl.MOCK[order(metaTbl.MOCK$rsum,decreasing=T),]
N=round(0.05*nrow(metaTbl.MOCK))
metaTbl.MOCK <- metaTbl.MOCK[-c(1:N,(nrow(metaTbl.MOCK)-N):nrow(metaTbl.MOCK)), ]
# metaTbl.MOCK[c(1:120)] <- metaTbl.MOCK[c(1:120)]/metaTbl.MOCK$rsum
metaTbl.MOCK$rsum <- NULL

metaTbl.JTE$rsum <- apply(metaTbl.JTE, 1, sum,na.rm = T)
metaTbl.JTE <- metaTbl.JTE[order(metaTbl.JTE$rsum,decreasing=T),]
N=round(0.05*nrow(metaTbl.JTE))
metaTbl.JTE <- metaTbl.JTE[-c(1:N,(nrow(metaTbl.JTE)-N):nrow(metaTbl.JTE)), ]
# metaTbl.JTE[c(1:120)] <- metaTbl.JTE[c(1:120)]/metaTbl.JTE$rsum
metaTbl.JTE$rsum <- NULL






plotTbl <- data.frame("DMSO"=apply(metaTbl.MOCK, 2, mean, na.rm = T),"JTE"=apply(metaTbl.JTE, 2, mean, na.rm = T))
plotTbl$pos <- c(1:120)
# plotTbl$DMSO <- as.numeric(CMA(plotTbl$DMSO))
# plotTbl$JTE <- as.numeric(CMA(plotTbl$JTE))
plotTbl <- melt(plotTbl,id.vars = "pos")
names(plotTbl) <- c("pos", "sample", "mean.rpm")
# write.csv(plotTbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/TES_U937MP.csv"),row.names = F)

plot[[paste0("U937MP",".TES")]] <- ggplot(plotTbl, aes(x=pos, y=mean.rpm, group=sample, color=sample)) +
  geom_line(size=1)+
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  scale_x_continuous(limits = c(10,110),breaks=seq(10,110,50),labels=c("-5kb","TES","5kb")) +
  scale_color_manual(values=c("grey","red")) + 
  coord_cartesian(ylim = c(0.05, 0.22)) +
  scale_y_continuous(breaks = seq(0.05, 0.2, by=0.025)) +
  labs(title = "U937MP")





















library(reshape)
library(tidyverse)
`%notin%` <- Negate(`%in%`)


for (cell in c("U937","HepG2","HeLa")) {
  
  # cell="U937"
  
  
  metaTbl.plus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-08/_metagene/",cell,"_maxSize.plus.tab"), header = T, sep = "\t")
  metaTbl.minus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-08/_metagene/",cell,"_maxSize.minus.tab"), header = T, sep = "\t")
  names(metaTbl.plus) <- gsub(".plus","",names(metaTbl.plus))
  names(metaTbl.minus) <- gsub(".minus","",names(metaTbl.minus))
  
  metaTbl <- rbind(metaTbl.plus,metaTbl.minus)
  
  metaTbl.MOCK <- metaTbl[,grep("DMSO",names(metaTbl))]
  metaTbl.JTE <- metaTbl[,grep("JTE",names(metaTbl))]
  
  
  metaTbl.MOCK$rsum <- apply(metaTbl.MOCK, 1, sum,na.rm = T)
  metaTbl.MOCK <- metaTbl.MOCK[order(metaTbl.MOCK$rsum,decreasing=T),]
  N=round(0.05*nrow(metaTbl.MOCK))
  metaTbl.MOCK <- metaTbl.MOCK[-c(1:N,(nrow(metaTbl.MOCK)-N):nrow(metaTbl.MOCK)), ]
  # metaTbl.MOCK[c(1:220)] <- metaTbl.MOCK[c(1:220)]/metaTbl.MOCK$rsum
  metaTbl.MOCK$rsum <- NULL
  
  metaTbl.JTE$rsum <- apply(metaTbl.JTE, 1, sum,na.rm = T)
  metaTbl.JTE <- metaTbl.JTE[order(metaTbl.JTE$rsum,decreasing=T),]
  N=round(0.05*nrow(metaTbl.JTE))
  metaTbl.JTE <- metaTbl.JTE[-c(1:N,(nrow(metaTbl.JTE)-N):nrow(metaTbl.JTE)), ]
  # metaTbl.JTE[c(1:220)] <- metaTbl.JTE[c(1:220)]/metaTbl.JTE$rsum
  metaTbl.JTE$rsum <- NULL
  
  
  
  
  
  
  plotTbl <- data.frame("DMSO"=apply(metaTbl.MOCK, 2, mean, na.rm = T),"JTE"=apply(metaTbl.JTE, 2, mean, na.rm = T))
  plotTbl$pos <- c(1:220)
  # plotTbl$DMSO <- as.numeric(CMA(plotTbl$DMSO))
  # plotTbl$JTE <- as.numeric(CMA(plotTbl$JTE))
  plotTbl <- melt(plotTbl,id.vars = "pos")
  names(plotTbl) <- c("pos", "sample", "mean.rpm")
  # write.csv(plotTbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/whole_",cell,".csv"),row.names = F)
  
  plot[[paste0(cell,".whole")]] <- ggplot(plotTbl, aes(x=pos, y=mean.rpm, group=sample, color=sample)) +
    geom_line(size=1)+
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    scale_x_continuous(limits = c(10,210),breaks=seq(10,210,50),labels=c("-5kb","TSS","","TES","5kb")) +
    scale_color_manual(values=c("grey","red")) + 
    coord_cartesian(ylim = c(0.07, 0.17)) +
    scale_y_continuous(breaks = seq(0.07, 0.17, by=0.02)) +
    labs(title = cell)
  
  
  
  
  
  
  
  metaTbl.plus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-08/_metagene/TES5K.",cell,"_maxSize.plus.tab"), header = T, sep = "\t")
  metaTbl.minus <- read.table(paste0(workplace,"/linux/lwang/drive/project/21029-08/_metagene/TES5K.",cell,"_maxSize.minus.tab"), header = T, sep = "\t")
  names(metaTbl.plus) <- gsub(".plus","",names(metaTbl.plus))
  names(metaTbl.minus) <- gsub(".minus","",names(metaTbl.minus))
  
  metaTbl <- rbind(metaTbl.plus,metaTbl.minus)
  
  metaTbl.MOCK <- metaTbl[,grep("DMSO",names(metaTbl))]
  metaTbl.JTE <- metaTbl[,grep("JTE",names(metaTbl))]
  
  
  metaTbl.MOCK$rsum <- apply(metaTbl.MOCK, 1, sum,na.rm = T)
  metaTbl.MOCK <- metaTbl.MOCK[order(metaTbl.MOCK$rsum,decreasing=T),]
  N=round(0.05*nrow(metaTbl.MOCK))
  metaTbl.MOCK <- metaTbl.MOCK[-c(1:N,(nrow(metaTbl.MOCK)-N):nrow(metaTbl.MOCK)), ]
  # metaTbl.MOCK[c(1:120)] <- metaTbl.MOCK[c(1:120)]/metaTbl.MOCK$rsum
  metaTbl.MOCK$rsum <- NULL
  
  metaTbl.JTE$rsum <- apply(metaTbl.JTE, 1, sum,na.rm = T)
  metaTbl.JTE <- metaTbl.JTE[order(metaTbl.JTE$rsum,decreasing=T),]
  N=round(0.05*nrow(metaTbl.JTE))
  metaTbl.JTE <- metaTbl.JTE[-c(1:N,(nrow(metaTbl.JTE)-N):nrow(metaTbl.JTE)), ]
  # metaTbl.JTE[c(1:120)] <- metaTbl.JTE[c(1:120)]/metaTbl.JTE$rsum
  metaTbl.JTE$rsum <- NULL
  
  
  
  
  
  
  plotTbl <- data.frame("DMSO"=apply(metaTbl.MOCK, 2, mean, na.rm = T),"JTE"=apply(metaTbl.JTE, 2, mean, na.rm = T))
  plotTbl$pos <- c(1:120)
  # plotTbl$DMSO <- as.numeric(CMA(plotTbl$DMSO))
  # plotTbl$JTE <- as.numeric(CMA(plotTbl$JTE))
  plotTbl <- melt(plotTbl,id.vars = "pos")
  names(plotTbl) <- c("pos", "sample", "mean.rpm")
  # write.csv(plotTbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/TES_",cell,".csv"),row.names = F)
  
  plot[[paste0(cell,".TES")]] <- ggplot(plotTbl, aes(x=pos, y=mean.rpm, group=sample, color=sample)) +
    geom_line(size=1)+
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    scale_x_continuous(limits = c(10,110),breaks=seq(10,110,50),labels=c("-5kb","TES","5kb")) +
    scale_color_manual(values=c("grey","red")) + 
    coord_cartesian(ylim = c(0.05, 0.22)) +
    scale_y_continuous(breaks = seq(0.05, 0.2, by=0.025)) +
    labs(title = cell)

  
}


library(cowplot)
plot_grid(plotlist=plot,ncol = 2)
plot_grid(plotlist=plot[c(3)],ncol = 1)
plot_grid(plotlist=plot[c(5)],ncol = 1)

# plot_grid(plotlist=plot[c(2,4,6,8)],ncol = 2)

