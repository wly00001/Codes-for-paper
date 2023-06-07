#####################
### input DE data ###

library(reshape)
library(tidyverse)
library(ggrepel)
library(cowplot)

# mychart.Correlation <- function (R, method = c("pearson", "kendall", "spearman"), ...){
#   
#   x = PerformanceAnalytics::checkData(R, method = "matrix")
#   if (missing(method)) method = method[1]
#   cormeth <- method
#   panel.cor <- function(x, y, digits = 2, prefix = "", 
#                         use = "pairwise.complete.obs", 
#                         method = cormeth, 
#                         cex.cor, ...) {
#     usr <- par("usr")
#     on.exit(par(usr))
#     par(usr = c(0, 1, 0, 1))
#     r <- cor(x, y, use = use, method = method)
#     txt <- format(c(r, 0.123456789), digits = digits)[1]
#     txt <- paste(prefix, txt, sep = "")
#     if (missing(cex.cor)) cex <- 1.5
#     test <- cor.test(as.numeric(x), as.numeric(y), method = method)
#     Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
#                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
#                      symbols = c("***", "**", "*", ".", " "))
#     text(0.5, 0.5, txt, cex = cex)
#     text(0.8, 0.8, Signif, cex = cex, col = 2)
#   }
#   
#   pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
#   
# }




##########################
### gene feature table ###

genelist <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/gene_features_hg19/hg19_refseq_gene_features_20220707.txt", header = T, quote = "", sep = "\t")[c(1,4,2,3,9,12)]
genelist <- unique(genelist)
RNA_stability <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/gene_features_hg19/RNA_stability_2022Aysegul/RNA_stability_2022_Aysegul.tbl", header = T, quote = "", sep = "\t")[c(1,3)]
names(RNA_stability) <- c("gene_symbol","RNA_stability_HepG2")

genelist <- merge(genelist, RNA_stability)


gene_density_100kb <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/gene_features_hg19/hg19_refseq_gene_density_100kb.txt", header = T, quote = "", sep = "\t")
names(gene_density_100kb) <- c("gene_symbol","gene_density_100kb")
# gene_density_1mb <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/gene_features_hg19/hg19_refseq_gene_density_1mb.txt", header = T, quote = "", sep = "\t")
# names(gene_density_1mb) <- c("gene_symbol","gene_density_1mb")
# gene_density_10mb <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/gene_features_hg19/hg19_refseq_gene_density_10mb.txt", header = T, quote = "", sep = "\t")
# names(gene_density_10mb) <- c("gene_symbol","gene_density_10mb")


genelist <- merge(genelist, gene_density_100kb)
# genelist <- merge(genelist, gene_density_1mb)
# genelist <- merge(genelist, gene_density_10mb)

genelist[c(3:6,8)]= log10(genelist[c(3:6,8)])
names(genelist)[c(3:6,8)] <- paste0("log10.",names(genelist)[c(3:6,8)])
genelist <- genelist %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))  
genelist <- na.omit(genelist)
genelist <- genelist[!duplicated(genelist$gene_symbol),]



tbl <- data.frame(genelist[-c(1)],row.names=genelist$gene_symbol)

# mychart.Correlation(tbl, histogram=FALSE, method = "pearson")


# cor_all <- cor.test(tbl$gene_GC,tbl$RNA_stability_HepG2)
# ggplot(data=tbl, aes(x=gene_GC, y=RNA_stability_HepG2)) +
#   geom_point() +
#   geom_hex(bins = 75) +
#   # coord_fixed(ratio = 1)+
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         # legend.position = "none",
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   ggtitle(paste0(#"REDi, ",sample," vs neuron",
#                  "\ncor=",round(cor_all$estimate,2),", pval=",formatC(cor_all$p.value, format = "e", digits = 2)))
# 
# cor_all <- cor.test(tbl$gene_GC,tbl$gene_density_100kb)
# ggplot(data=tbl, aes(x=gene_GC, y=gene_density_100kb)) +
#   geom_point() +
#   geom_hex(bins = 75) +
#   # coord_fixed(ratio = 1)+
#   theme_classic()+
#   theme(axis.text.y=element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text.x=element_text(colour = "black"),
#         axis.title.x = element_text(colour = "black"),
#         # legend.position = "none",
#         axis.line = element_line(colour = "black"),
#         axis.ticks = element_line(size = 1),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
#   ggtitle(paste0(#"REDi, ",sample," vs neuron",
#     "\ncor=",round(cor_all$estimate,2),", pval=",formatC(cor_all$p.value, format = "e", digits = 2)))


library(tidyverse)
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

bin_tbl <- data.frame("Bin_no."=c("bin1","bin2","bin3","bin4","bin5"))

DF <- tbl[order(tbl$gene_GC,decreasing=F),]
bin1=bin2=bin3=bin4=round(nrow(DF)/5)
bin5=nrow(DF)-(bin1+bin2+bin3+bin4)
DF$bin <- c(rep("bin1",bin1),rep("bin2",bin2),rep("bin3",bin3),rep("bin4",bin4),rep("bin5",bin5))

bin_tbl2<-data.frame("gene_N"=c(nrow(subset(DF, bin == "bin1")),
                                    nrow(subset(DF, bin == "bin2")),
                                    nrow(subset(DF, bin == "bin3")),
                                    nrow(subset(DF, bin == "bin4")),
                                    nrow(subset(DF, bin == "bin5"))),
                         "GC_range"=c(paste(range(subset(DF, bin == "bin1")$gene_GC),collapse = "-"),
                                            paste(range(subset(DF, bin == "bin2")$gene_GC),collapse = "-"),
                                            paste(range(subset(DF, bin == "bin3")$gene_GC),collapse = "-"),
                                            paste(range(subset(DF, bin == "bin4")$gene_GC),collapse = "-"),
                                            paste(range(subset(DF, bin == "bin5")$gene_GC),collapse = "-")))
bin_tbl <- cbind(bin_tbl,bin_tbl2)


GC_bin_table <- list()
boxplot <- list()
for (c in c("log10.gene_size","log10.max_tx_size","log10.overall_intron_size","log10.max_intron_size","RNA_stability_HepG2","log10.gene_density_100kb")){
  
  # c="log10.gene_size"
  
  DF2 <- DF[c(c,"bin")]
  names(DF2) <- c("feature","bin")
  
  DF3 <- data_summary(DF2, varname="feature", 
                      groupnames=c("bin"))
  DF3$featureName <- rep(c, nrow(DF3))
  GC_bin_table[[c]] <- DF3
  
  KS <- ks.test(subset(DF2, bin == "bin1")$feature,subset(DF2, bin == "bin5")$feature)
  # ggplot(DF3, aes(x=bin, y=feature, group=1)) + 
  #   geom_line(size=1) +
  #   geom_point(size=3) +
  #   geom_errorbar(aes(ymin=feature-se, ymax=feature+se), width=0.5, size=1,
  #                 position=position_dodge(0.05)) +
  #   theme_classic()+
  #   theme(axis.text.y=element_text(colour = "black"),
  #         axis.title.y = element_text(colour = "black"),
  #         axis.text.x=element_text(colour = "black"),
  #         axis.title.x = element_text(colour = "black"),
  #         axis.line = element_line(colour = "black"),
  #         axis.ticks = element_line(size = 1),
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  #   labs(title = paste0("KS_bin1.v.bin5=",round(KS$statistic,2),", pval=",formatC(KS$p.value, format = "e", digits = 2)),
  #        x = "GC content bin", y = "Mean of RNA_stability_HepG2")
  
  
  boxplot[[c]] <- ggplot(DF2, aes(x=bin, y=feature,fill=bin)) + 
    geom_boxplot(outlier.shape = NA,coef=0) + 
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
    coord_cartesian(ylim = c(min(quantile(subset(DF2,bin=="bin1")$feature, 0.1),
                                      quantile(subset(DF2,bin=="bin2")$feature, 0.1),
                                      quantile(subset(DF2,bin=="bin3")$feature, 0.1),
                                      quantile(subset(DF2,bin=="bin4")$feature, 0.1),
                                      quantile(subset(DF2,bin=="bin5")$feature, 0.1)),
                                  max(quantile(subset(DF2,bin=="bin1")$feature, 0.9),
                                      quantile(subset(DF2,bin=="bin2")$feature, 0.9),
                                      quantile(subset(DF2,bin=="bin3")$feature, 0.9),
                                      quantile(subset(DF2,bin=="bin4")$feature, 0.9),
                                      quantile(subset(DF2,bin=="bin5")$feature, 0.9)))) +
    labs(title = paste0("KS_bin1.v.bin5=",round(KS$statistic,2),", pval=",formatC(KS$p.value, format = "e", digits = 2)),
         x = "GC content bin", y = c)
  
  
  
}

library(cowplot)
plot_grid(plotlist=boxplot,ncol = 2)


GC_bin_table2 <- do.call(rbind,GC_bin_table)
# write.csv(GC_bin_table2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig3c.csv"),row.names = F)


