# REDu measures the relative expression levels between the top two most abundant APA isoforms in the 3'-most exon only.
# REDi measures the relative expression levels between the top expressed isoform in the 3'-most exon and the top expressed isoform in an intron.



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




library(reshape2)
library(tidyverse)
library(ggrepel)

# PAS_aaa1 <- PAS_aaa <- read.csv(paste0("./pas_U937.csv"), header = T)

List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_HeLa_HepG2/List.REDu.RData"))[c(1:4)]
names(List.REDu) <- c("HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h")

#REDi


List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x<-na.omit(x); x })




List.REDu2 <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDu2[[sample]] <- df
  
}

List.REDu.8h3 <- lapply(List.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu2 <- heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3)
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu2, heatmapTbl.REDu2[[3]]!="NO" | heatmapTbl.REDu2[[5]]!="NO" | heatmapTbl.REDu2[[7]]!="NO" | heatmapTbl.REDu2[[9]]!="NO")



################################################################################
################################################################################
################################################################################

# aUTR_size_vs_REDu

dir.create(file.path("./", "_aUTR_size_vs_REDu"))


plotTbl <-list()
bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

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

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  REDu_aaa <- List.REDu[[sample]]
  REDu_aaa <- subset(REDu_aaa, gene_symbol %in% heatmapTbl.REDu2$gene_symbol)
  
  is.na(REDu_aaa) <- sapply(REDu_aaa, is.infinite)
  REDu_aaa <- na.omit(REDu_aaa)
  REDu_aaa <- REDu_aaa[order(REDu_aaa$aUTR_size,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(REDu_aaa)/5)
  bin5=nrow(REDu_aaa)-(bin1+bin2+bin3+bin4)
  REDu_aaa$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  
  bin_tbl2<-data.frame("gene_N"=c(nrow(subset(REDu_aaa, bin == "bin1")),
                                  nrow(subset(REDu_aaa, bin == "bin2")),
                                  nrow(subset(REDu_aaa, bin == "bin3")),
                                  nrow(subset(REDu_aaa, bin == "bin4")),
                                  nrow(subset(REDu_aaa, bin == "bin5"))),
                       "aUTR_size_range"=c(paste(range(subset(REDu_aaa, bin == "bin1")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin2")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin3")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin4")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin5")$aUTR_size),collapse = "~")))
  names(bin_tbl2) <- paste(names(bin_tbl2),sample,sep = "_")
  bin_tbl <- cbind(bin_tbl,bin_tbl2)
  
  REDu_aaa3 <- data_summary(REDu_aaa, varname="REDu1", 
                            groupnames=c("bin"))
  REDu_aaa3$sampleName <- rep(sample, nrow(REDu_aaa3))
  plotTbl[[sample]] <- REDu_aaa3
  
  
  
}

plotTbl_JTE <- do.call(rbind,plotTbl)

library(RColorBrewer)
ggplot(plotTbl_JTE, aes(x=bin, y=REDu1, group=sampleName, color=sampleName)) + 
  geom_line(size=1) +
  geom_point(size=3) +
  # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
  geom_errorbar(aes(ymin=REDu1-se, ymax=REDu1+se), width=0.5, size=1,
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
  labs(x = "aUTR size bin", y = "Mean of REDu") +
  scale_color_manual(values=c(colorHeLa_1,colorHeLa_10,colorHepG2_1,colorHepG2_10))
# write.csv(as.data.frame(t(bin_tbl)), "./_aUTR_size_vs_REDu/bin_tbl.csv", col.names = F, row.names = T)
# write.csv(plotTbl_JTE, "./_aUTR_size_vs_REDu/bin_plot_tbl.csv", col.names = F, row.names = T)

# write.csv(plotTbl_JTE, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig6a.csv"),row.names = T)










# REDu measures the relative expression levels between the top two most abundant APA isoforms in the 3'-most exon only.
# REDi measures the relative expression levels between the top expressed isoform in the 3'-most exon and the top expressed isoform in an intron.

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

library(reshape2)
library(tidyverse)
library(ggrepel)

List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum1011paper/APA_U937_U937P/List.REDu.RData"))[c(1:2)]
names(List.REDu) <- c("U937_1_7h","U937PMA_1_7h")

#REDi


List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x<-na.omit(x); x })




List.REDu2 <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDu2[[sample]] <- df
  
}

List.REDu.8h3 <- lapply(List.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3)
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO")



################################################################################
################################################################################
################################################################################

# aUTR_size_vs_REDu

dir.create(file.path("./", "_aUTR_size_vs_REDu"))


plotTbl <-list()
bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

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

dataTbl <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  REDu_aaa <- List.REDu[[sample]]
  REDu_aaa <- subset(REDu_aaa, gene_symbol %in% heatmapTbl.REDu2$gene_symbol)
  
  is.na(REDu_aaa) <- sapply(REDu_aaa, is.infinite)
  REDu_aaa <- na.omit(REDu_aaa)
  REDu_aaa <- REDu_aaa[order(REDu_aaa$aUTR_size,decreasing=F),]
  bin1=bin2=bin3=bin4=round(nrow(REDu_aaa)/5)
  bin5=nrow(REDu_aaa)-(bin1+bin2+bin3+bin4)
  REDu_aaa$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))
  dataTbl[[sample]] <- REDu_aaa
  
  bin_tbl2<-data.frame("gene_N"=c(nrow(subset(REDu_aaa, bin == "bin1")),
                                  nrow(subset(REDu_aaa, bin == "bin2")),
                                  nrow(subset(REDu_aaa, bin == "bin3")),
                                  nrow(subset(REDu_aaa, bin == "bin4")),
                                  nrow(subset(REDu_aaa, bin == "bin5"))),
                       "aUTR_size_range"=c(paste(range(subset(REDu_aaa, bin == "bin1")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin2")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin3")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin4")$aUTR_size),collapse = "~"),
                                           paste(range(subset(REDu_aaa, bin == "bin5")$aUTR_size),collapse = "~")))
  names(bin_tbl2) <- paste(names(bin_tbl2),sample,sep = "_")
  bin_tbl <- cbind(bin_tbl,bin_tbl2)
  
  REDu_aaa3 <- data_summary(REDu_aaa, varname="REDu1", 
                            groupnames=c("bin"))
  REDu_aaa3$sampleName <- rep(sample, nrow(REDu_aaa3))
  plotTbl[[sample]] <- REDu_aaa3
  
  
  
}

plotTbl_JTE <- do.call(rbind,plotTbl)

library(RColorBrewer)
ggplot(plotTbl_JTE, aes(x=bin, y=REDu1, group=sampleName, color=sampleName)) + 
  geom_line(size=1) +
  geom_point(size=3) +
  # scale_color_manual(labels = c("U937_JTE607_6h/U937_Mock_6h", "U937_PMA_JTE607_6h/U937_PMA_Mock_6h"), values=c("grey","black")) +
  geom_errorbar(aes(ymin=REDu1-se, ymax=REDu1+se), width=0.5, size=1,
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
  labs(title = paste0("bin5 comparison, wilcox.test, ", formatC(wilcox.test(dataTbl[[1]]$REDu1, dataTbl[[2]]$REDu1)$p.value, format = "e", digits = 2)),
       x = "aUTR size bin", y = "Mean of REDu") +
  scale_color_manual(values=c(colorU937,colorU937MP))

# write.csv(as.data.frame(t(bin_tbl)), "./bin_tbl.csv", col.names = F, row.names = T)
# write.csv(plotTbl_JTE, "./bin_plot_tbl.csv", col.names = F, row.names = T)

# write.csv(plotTbl_JTE, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS6d.csv"),row.names = T)


