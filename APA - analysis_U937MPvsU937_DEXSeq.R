workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"


# REDu measures the relative expression levels between the top two most abundant APA isoforms in the 3'-most exon only.
# REDi measures the relative expression levels between the top expressed isoform in the 3'-most exon and the top expressed isoform in an intron.

library(reshape2)
library(tidyverse)
library(ggrepel)


filenames <- list.files(paste0(workplace,"/linux/lwang/JTEonly_APA_strPAS/U937_1h6h/result/U937_1h6h/DEXSeq_out_0310/"), pattern="*.tbl", full.names=TRUE)
filenames

REDulist <- list()
plot_REDu <- list()
ratio_tbl_REDu <- data.frame(comparison=character(),
                             TPA_LE=numeric(),
                             TPA_SH=numeric(),
                             TPA_Ratio=numeric(),
                             stringsAsFactors=FALSE)

for (file in filenames) {

  sample <- gsub(".tbl","",tail(strsplit(file,"/")[[1]],1))
  df <- read.table(paste0(file), header = T, sep = "\t")
  
  df$change <- "NO"
  df$change[df$pval_gene < 0.05 &  df$RED < (-log2(1.2))] <- "shortened"
  df$change[df$pval_gene < 0.05 &  df$RED > log2(1.2)] <- "lengthened"
  
  REDulist[[sample]] <- df
  
  plot_REDu[[sample]] <- ggplot() +
    geom_point(data=subset(df, change == "NO"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="NO"), size = 0.1, alpha=0.3)+
    geom_point(data=subset(df, change == "lengthened"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="lengthened"), size = 1.5, alpha=0.3)+
    geom_point(data=subset(df, change == "shortened"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="shortened"), size = 1.5, alpha=0.3)+
    scale_color_manual(values=c("red", "grey", "blue")) +
    theme_classic()+
    theme(axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text.x=element_text(colour = "black"),
          axis.title.x = element_text(colour = "black"),
          plot.title = element_text(size = 10),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1)) +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    scale_y_continuous(limits = c(-3, 3), breaks = c(-4,-2,0,2,4))+
    # annotate('text', x=-2, y=2, label=nrow(subset(df, change == "lengthened")),color = "red") +
    # annotate('text', x=2, y=-2, label=nrow(subset(df, change == "shortened")),color = "blue") +
    # geom_label_repel(data=subset(df, gene %in% c("TIMP2")),
    #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
    labs(title = paste(sample,
                       "\nred=",table(df$change)[1],
                       ", gray=",table(df$change)[2],
                       ", blue=",table(df$change)[3],
                       ", ratio= ",ifelse(table(df$change)[1]>table(df$change)[3],round(table(df$change)[1]/table(df$change)[3],1),(-1)*round(table(df$change)[3]/table(df$change)[1],1)),
                       sep = ""),
         x = "Log2Ratio of pPAS isoforms", y = "Log2Ratio of dPAS isoforms", color = "")

  print(table(df$change))
  ratio_tbl_REDu[nrow(ratio_tbl_REDu) + 1,] = c(sample,table(df$change)[1],table(df$change)[3],ifelse(table(df$change)[1]>table(df$change)[3],round(table(df$change)[1]/table(df$change)[3],1),(-1)*round(table(df$change)[3]/table(df$change)[1],1)))
  
}

library(cowplot)
plot_grid(plotlist=c(plot_REDu),ncol = 1)

# write.table(ratio_tbl_REDu, "REDu.tbl", row.names = F, col.names = T, sep = "\t", quote = F)
# write.csv(REDulist[[1]], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig7b.csv"),row.names = T)


# REDilist <- list()
# plot_REDi <- list()
# ratio_tbl_REDi <- data.frame(comparison=character(),
#                              IPA_SU=numeric(),
#                              IPA_AC=numeric(),
#                              IPA_Ratio=numeric(),
#                              stringsAsFactors=FALSE)
# 
# for (file in filenames[c(1:4)]) {
#   
#   sample <- gsub(".tbl","",tail(strsplit(file,"/")[[1]],1))
#   df <- read.table(paste0(file), header = T, sep = "\t")
#   
#   df <- df %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
#   df <- na.omit(df)
#   
#   df$change <- "NO"
#   df$change[df$pval_gene < 0.05 &  df$RED < (-log2(1.2))] <- "shortened"
#   df$change[df$pval_gene < 0.05 &  df$RED > log2(1.2)] <- "lengthened"
#   
#   REDilist[[sample]] <- df
#   
#   plot_REDi[[sample]] <- ggplot() +
#     geom_point(data=subset(df, change == "NO"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="NO"), size = 0.1, alpha=0.3)+
#     geom_point(data=subset(df, change == "lengthened"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="lengthened"), size = 1.5, alpha=0.3)+
#     geom_point(data=subset(df, change == "shortened"), aes(x=as.numeric(pPASratio), y=as.numeric(dPASratio), color="shortened"), size = 1.5, alpha=0.3)+
#     scale_color_manual(values=c("red", "grey", "blue")) +
#     theme_classic()+
#     theme(axis.text.y=element_text(colour = "black"),
#           axis.title.y = element_text(colour = "black"),
#           axis.text.x=element_text(colour = "black"),
#           axis.title.x = element_text(colour = "black"),
#           plot.title = element_text(size = 10),
#           legend.position = "none",
#           axis.line = element_line(colour = "black"),
#           axis.ticks = element_line(size = 1),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#     coord_fixed(ratio = 1)+
#     scale_x_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
#     scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
#     # annotate('text', x=-2, y=2, label=nrow(subset(df, change == "lengthened")),color = "red") +
#     # annotate('text', x=2, y=-2, label=nrow(subset(df, change == "shortened")),color = "blue") +
#     # geom_label_repel(data=subset(df, gene %in% c("TIMP2")),
#     #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
#     labs(title = paste(sample,
#                        "\nred=",table(df$change)[1],
#                        ", gray=",table(df$change)[2],
#                        ", blue=",table(df$change)[3],
#                        ", ratio= ",ifelse(table(df$change)[1]>table(df$change)[3],round(table(df$change)[1]/table(df$change)[3],1),(-1)*round(table(df$change)[3]/table(df$change)[1],1)),
#                        sep = ""),
#          x = "Log2Ratio of IPA isoforms", y = "Log2Ratio of TPA isoforms", color = "")
#   
#   print(table(df$change))
#   ratio_tbl_REDi[nrow(ratio_tbl_REDi) + 1,] = c(sample,table(df$change)[1],table(df$change)[3],ifelse(table(df$change)[1]>table(df$change)[3],round(table(df$change)[1]/table(df$change)[3],1),(-1)*round(table(df$change)[3]/table(df$change)[1],1)))
#   
# }
# 
# library(cowplot)
# plot_grid(plotlist=c(plot_REDi),ncol = 4)
# 
# # write.table(ratio_tbl_REDi, "REDi.tbl", row.names = F, col.names = T, sep = "\t", quote = F)
# 
# ratio_tbl3 <- cbind(ratio_tbl_REDu,ratio_tbl_REDi)
# 
# 
# write.table(ratio_tbl3, paste0("./APA ratio 0310.tbl"), row.names = F, col.names = T, sep = "\t")
# 
# 
# 
# 
# 
# 
# dir.create(file.path("./", "_APA_report_0310"))
# 
# NCBI = read.table("D:/wly/OneDrive - The Wistar Institute/project/reference/Homo_sapiens.gene_info", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = c("","NA"), fill = TRUE)
# NCBI = as.data.frame(sapply(NCBI, function(x) gsub("\"", "", x)))
# 
# 
# 
# 
# List.REDu <- lapply(REDulist, function(x) { x<-x[c(1,2,3,5,6,7,28,29)]; x })
# List.REDu <- lapply(List.REDu, function(x) { x<-merge(NCBI[c("Symbol","description")], x, all.y = T, by.x = "Symbol", by.y = "gene_symbol"); x })
# List.REDu <- lapply(List.REDu, function(x) { x[c(5,6)]<-format(x[c(5,6)],big.mark=",",scientific=F); x })
# List.REDu <- lapply(List.REDu, function(x) { x[c(8)]<-round(x[c(8)],2); x })
# List.REDu <- lapply(List.REDu, function(x) { x[c(7)] <- formatC(unlist(x[c(7)]), format = "e", digits = 2); x })
# 
# 
# List.REDi <- lapply(REDilist, function(x) { x<-x[c(1,2,23,24)]; x })
# List.REDi <- lapply(List.REDi, function(x) { x<-merge(NCBI[c("Symbol","description")], x, all.y = T, by.x = "Symbol", by.y = "gene_symbol"); x })
# List.REDi <- lapply(List.REDi, function(x) { x[c(4)]<-round(x[c(4)],2); x })
# List.REDi <- lapply(List.REDi, function(x) { x[c(3)] <- formatC(unlist(x[c(3)]), format = "e", digits = 2); x })
# 
# 
# 
# 
# 
# 
# library(openxlsx)
# wb = createWorkbook()
# 
# for (n in 1:length(List.REDu)) {
#   # n=1
#   name <- names(List.REDu)[n]
#   
#   # Write to Excel
#   addWorksheet(wb, name)
#   writeData(wb, sheet = n, List.REDu[[n]], rowNames = F, colNames = T, keepNA = T)
#   
# }
# 
# saveWorkbook(wb, "./_APA_report_0310/APA_report_REDu.xlsx", overwrite = T)
# 
# 
# 
# 
# library(openxlsx)
# wb = createWorkbook()
# 
# for (n in 1:length(List.REDi)) {
#   # n=1
#   name <- names(List.REDi)[n]
#   
#   # Write to Excel
#   addWorksheet(wb, name)
#   writeData(wb, sheet = n, List.REDi[[n]], rowNames = F, colNames = T, keepNA = T)
#   
# }
# 
# saveWorkbook(wb, "./_APA_report_0310/APA_report_REDi.xlsx", overwrite = T)






