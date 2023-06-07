
library(tidyverse)
library(RColorBrewer)
library(reshape2)

`%!in%` <- Negate(`%in%`)
hg2mm <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/HomoloGene/homologene_hg2mm.tbl", header = T, sep = "\t", quote = "")[c(4,9)]


DoG <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/reference/pnas.1711120114.sd01.txt", header = F, sep = "\t")[c(1)]
names(DoG) <- c("gene_symbol")

DoG2 <- merge(DoG,hg2mm,by.x="gene_symbol",by.y="mm.Gene_Symbol")



List.RS <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_RS_table/List.RS_UR_GB_DR_4K.RData")
List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF.RData")



DGE_CDF <- list()
library(openxlsx)
wb = createWorkbook()

for (cell in names(List.DGE_QSF)) {
  
  # cell="HepG2"
  
  DGE <- List.DGE_QSF[[cell]][c(1,4)]
  names(DGE) <- c("gene_symbol","DMSO_RPM")
  RS <- List.RS[[cell]][c(1,13)]
  
  DF <- merge(RS,DGE)
  DF$class <- ifelse(DF$gene_symbol %in% DoG2$hg.Gene_Symbol, "DoG", "other")
  
  DF <- DF[order(DF$DMSO_RPM,decreasing=F),]
  DF_DoG <- subset(DF,class=="DoG")

  n=round(nrow(DF)*0.01)
  List_other <- list()
  Ntbl <- data.frame(cell=character(),
                     seed=numeric(),
                     N_DoG=numeric(),
                     N_other=numeric(),
                     stringsAsFactors=FALSE)
  
  for (s in c(1:20)) {
    
    DF_other <- DF[FALSE,]
    
    for (i in seq(1,nrow(DF),n)) {
      
      # i=1
      
      if (i+n-1 >= nrow(DF)) {
        df1 <- DF[c(i:nrow(DF)),]
      } else {
        df1 <- DF[c(i:(i+n-1)),]
      }
      
      N_DoG <- nrow(subset(df1,class=="DoG"))
      N_other <- nrow(subset(df1,class=="other"))
      Ntbl[nrow(Ntbl) + 1,] = c(cell,s,N_DoG,N_other)
      
      df_other <- subset(df1,class!="DoG")
      
      if (N_DoG > nrow(df1)/2) {
        set.seed(s)
        DF_other <- rbind(DF_other,df_other[sample(nrow(df_other), N_DoG, replace = T), ])
      } else {
        set.seed(s)
        DF_other <- rbind(DF_other,df_other[sample(nrow(df_other), N_DoG), ])
      }
      df2 <- DF_other[c(2)]
      df2 <- as.data.frame(df2[order(df2[[c(1)]],decreasing=F),])
      names(df2) <- paste0("deltaRTS_seed_",s)
      
      List_other[[s]] <- df2
      
    }
    
  }
  
  name <- cell
  
  # Write to Excel
  addWorksheet(wb, name)
  writeData(wb, sheet = match(cell,names(List.DGE_QSF)), Ntbl, rowNames = F, colNames = T)
  
  Tbl_other <- do.call(cbind,List_other)
  
  # write.csv(cbind(DF_DoG,Tbl_other), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS2b.csv"),row.names = F)
  
  
  Tbl_other$median <- apply(Tbl_other, 1, median)
  
  
  
  # DGE_CDF[[cell]] <- ggplot() +
  #   stat_ecdf(data=DF_DoG, aes(x=RS.delta.JTEvDMSO, color = "a"),geom = "step", size=1) +
  #   stat_ecdf(data=Tbl_other, aes(x=median, color = "b"),geom = "step", size=1) +
  #   scale_color_manual(values=c("black","grey"),
  #                      labels = c(paste0("DoG, ",nrow(DF_DoG),", ",round(median(DF_DoG$RS.delta.JTEvDMSO),2)),
  #                                 paste0("others, ",nrow(Tbl_other),", ",round(median(Tbl_other$median),2)))) +
  #   theme_classic()+
  #   theme(axis.text.y=element_text(colour = "black"),
  #         axis.title.y = element_text(colour = "black"),
  #         axis.text.x=element_text(colour = "black"),
  #         axis.title.x = element_text(colour = "black"),
  #         # legend.position = "none",
  #         plot.title = element_text(size = 10),
  #         axis.line = element_line(colour = "black"),
  #         axis.ticks = element_line(size = 1),
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  #   coord_cartesian(xlim = c(-1, 3.5)) +
  #   labs(title = paste0(cell, "\ndelta RTS, other genes vs DoG, DMSO RPM at similar level",
  #                       "\nDoGvOTHERS, ", formatC(ks.test(DF_DoG$RS.delta.JTEvDMSO, Tbl_other$median)$p.value, format = "e", digits = 2)),
  #        x = "delta RTS", y = "Cumulative fraction", color = "")
  
  
  
  DGE_CDF[[cell]] <- ggplot() +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_1, color = "a01"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_2, color = "a02"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_3, color = "a03"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_4, color = "a04"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_5, color = "a05"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_6, color = "a06"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_7, color = "a07"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_8, color = "a08"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_9, color = "a09"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_10, color = "a10"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_11, color = "a11"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_12, color = "a12"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_13, color = "a13"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_14, color = "a14"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_15, color = "a15"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_16, color = "a16"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_17, color = "a17"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_18, color = "a18"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_19, color = "a19"),geom = "step", size=1) +
    stat_ecdf(data=Tbl_other, aes(x=deltaRTS_seed_20, color = "a20"),geom = "step", size=1) +
    stat_ecdf(data=DF_DoG, aes(x=RS.delta.JTEvDMSO, color = "b"),geom = "step", size=1) +
    scale_color_manual(values=c(rep("grey",20),"black")) +
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          # legend.position = "none",
          plot.title = element_text(size = 10),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    coord_cartesian(xlim = c(-0.5, 3.5)) +
    labs(title = paste0(cell, "\ndelta RTS, other genes vs DoG, DMSO RPM at similar level",
                        "\nN=",nrow(DF_DoG),", ","DoG, ",round(median(DF_DoG$RS.delta.JTEvDMSO),2),", others, ",round(median(Tbl_other$median),2),
                        "\nDoGvOTHERS, ", formatC(ks.test(DF_DoG$RS.delta.JTEvDMSO, Tbl_other$median)$p.value, format = "e", digits = 2)),
         x = "delta RTS", y = "Cumulative fraction", color = "")
  
  
  
  
}

# saveWorkbook(wb, "./readthru - all vs panDoG - numbers.xlsx", overwrite = T)

library(cowplot)
plot_grid(plotlist=DGE_CDF,ncol = 2)
# plot_grid(plotlist=DGE_CDF[c(2)],ncol = 1) #for paper
