workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

library(DESeq2)
library(reshape)
library(tidyverse)
library(ggrepel)


gene_full_name <- read.table(paste0(workplace,"/linux/lwang/drive/project/18077-13/gene_full_names + additional info.txt"), header = T, sep = "\t", quote = "", na.strings = c("","NA"))
# gene_full_name <- gene_full_name[-c(1,2)]
gene_full_name <- gene_full_name %>% separate(gene_name, c("gene_name",NA),sep=",,")
gene_full_name <- gene_full_name[!duplicated(gene_full_name$gene_symbol), ]

id_convert <- read.table(paste0(workplace,"/linux/lwang/drive/project/18077-13/featureCounts/id convert.tbl"), header = F, sep = "\t")
names(id_convert) <- c("gene_id","gene_symbol")
id_convert2 <- id_convert[!duplicated(id_convert), ]

NCBI = read.table(paste0(workplace,"/linux/lwang/drive/project/reference/Homo_sapiens.gene_info"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, na.strings = c("","NA"), fill = TRUE)
NCBI = as.data.frame(sapply(NCBI, function(x) gsub("\"", "", x)))


DGElist <- list()

filenames <- list.files("./featureCounts_multi_str_QSF/", pattern="*.featurecounts.Rmatrix.txt", full.names=TRUE)[c(12,13,16,17)]
filenames
ldf <- lapply(filenames, read.table, header = T)
names(ldf) <- substr(filenames,31,nchar(filenames)-30)

ldf2 <- list()
for (n in 1:length(ldf)) {
  
  name <- names(ldf)[c(n)]
  
  df <- ldf[[n]]
  
  names(df) <- c("gene_id",name)
  
  ldf2[[name]] <- df
  
  
  
}


### OE

gene_expression <- Reduce(function(x, y) merge(x, y), c(list(id_convert2,gene_full_name),ldf2))
gene_expression_nonMT <- subset(gene_expression, Chromosome != "MT")
gene_expression2 <- gene_expression_nonMT[-c(1,3:8)]
gene_expression3 <- aggregate(. ~ gene_symbol, gene_expression2, sum)
gene_expression4 <- data.frame(gene_expression3[-c(1)],row.names=gene_expression3$gene_symbol)
gene_expression4[c(1:4)] <- round(gene_expression4[c(1:4)])
gene_expression5 <- gene_expression4 %>% filter_all(all_vars(. > 5))

colData <- data.frame("condition" = c(rep("U937",2),rep("U937MP",2)),
                      "type" = c(rep(NA,4)),
                      row.names = names(gene_expression5))
colData
dds1 <- dds <- DESeqDataSetFromMatrix(countData = gene_expression5, colData = colData, design = ~ condition)

dds1 <- estimateSizeFactors(dds1)
gene_expression6 <- data.frame(counts(dds1, normalized=TRUE))
names(gene_expression6) <- paste0("Norm_",names(gene_expression6))
gene_expression6$gene_symbol <- row.names(gene_expression6)

dds <- DESeq(dds)
dds


ctrl_samples=c("U937")
test_samples=c("U937MP")

resList <- list()
GEplot_MA <- list()
GEplot_volcano <- list()

for (k in 1:length(test_samples)) {
  
  # k=1
  
  
  ctrl <- ctrl_samples[k]
  test <- test_samples[k]
  
  res <- results(dds, contrast=c("condition",test,ctrl), alpha=0.05, independentFiltering=FALSE)
  resDF <- data.frame(res)
  resDF$gene_symbol <- row.names(resDF)
  resDF <- resDF[c(7,1,2,5,6)]
  resDF$regu.adj <- ifelse(resDF$padj < 0.05 & resDF$log2FoldChange>1, "UP", ifelse(resDF$padj < 0.05 & resDF$log2FoldChange<(-1), "DN", "NO"))
  resDF <- merge(NCBI[c("Symbol","description")], resDF, all.y = T, by.x = "Symbol", by.y = "gene_symbol")
  names(resDF)[c(1)] <- "gene_symbol"
  
  write.csv(resDF, paste0("./_DGE_RPM/QSF.",test,".v.",ctrl,".DESeq2.csv"), row.names = FALSE)
  
  resList[[paste0(test,".v.",ctrl)]] <- resDF
  
  GEplot_MA[[paste0(test,".v.",ctrl)]] <- ggplot(data = resDF, aes(x = log10(baseMean), y = log2FoldChange, col = regu.adj)) +
    geom_point(alpha=0.3) +
    scale_color_manual(values=c("blue", "grey", "red")) +
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
    # geom_label_repel(data=subset(resDF, gene_symbol %in% c("CG15160")),
    #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
    coord_cartesian(xlim = c(1, 5), ylim = c(-4, 4)) +
    ggtitle(paste("DGE, ",test,".v.",ctrl,
                  "\nDN=",table(resDF$regu.adj)[1],
                  ", NC=",table(resDF$regu.adj)[2],
                  ", UP=",table(resDF$regu.adj)[3],
                  ", UP/DN=",round(table(resDF$regu.adj)[3]/table(resDF$regu.adj)[1],2),sep = ""))
  
  GEplot_volcano[[paste0(test,".v.",ctrl)]] <- ggplot(data=resDF, aes(x=log2FoldChange, y=-log10(padj), col=regu.adj)) +
    geom_point(alpha=0.3) +
    scale_color_manual(values=c("blue", "grey", "red")) +
    geom_vline(xintercept=c(-log2(1.2), log2(1.2)), col="red", linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
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
    ggtitle(paste("DGE, ",test,".v.",ctrl,
                  "\nDN=",table(resDF$regu.adj)[1],
                  ", NC=",table(resDF$regu.adj)[2],
                  ", UP=",table(resDF$regu.adj)[3],
                  ", UP/DN=",round(table(resDF$regu.adj)[3]/table(resDF$regu.adj)[1],2),sep = ""))
  
  
}

library(cowplot)
plot_grid(plotlist=GEplot_MA,ncol = 1)
plot_grid(plotlist=GEplot_volcano,ncol = 1)

# write.csv(resDF, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS7b.csv"),row.names = T)









################################################################################
################################################################################
################################################################################

library(tidyverse)
library(RColorBrewer)

filenames <- list.files("./_DGE_RPM/FC2.0/", pattern="*.csv", full.names=TRUE)[c(1:12)]
filenames

GEplot_MA <- list()

for (k in 1:length(filenames)) {
  
  # k=1
  
  treat_gene_expression5 <- read.table(filenames[k], header = T, sep = ",")[c(1,4:6,9,10)]
  treat_gene_expression5$meanRPM <- (treat_gene_expression5[[2]]+treat_gene_expression5[[3]])/2
  treat_gene_expression5[c(2,3)] <-NULL
  names(treat_gene_expression5) <- c("gene_symbol","L2FC","pval.fisher.adj","regu.adj","meanRPM")
  
  # treat_gene_expression5[["regu.adj"]] <- ifelse(treat_gene_expression5[[paste0("pval.fisher.adj")]] < 0.05 & treat_gene_expression5[[paste0("L2FC")]]>log2(2), "UP", ifelse(treat_gene_expression5[[paste0("pval.fisher.adj")]] < 0.05 & treat_gene_expression5[[paste0("L2FC")]]<(-log2(2)), "DN", "NO")) #set log2FC cutoff as log2(2)
  
  sample <- gsub(".pseudocount_BH.csv","",tail(strsplit(filenames[k],"/")[[1]],1))
  
  ###volcano plots to show regulated genes
  GEplot_MA[[sample]] <- ggplot(data = treat_gene_expression5, aes(x = log2(meanRPM), y = L2FC, col = regu.adj)) +
    geom_point(alpha=0.3) +
    scale_color_manual(values=c("blue", "grey", "red")) +
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
    # geom_label_repel(data=subset(treat_gene_expression5, gene_symbol %in% c("CG15160")),
    #                  box.padding = 1.5,point.padding = 0.5,color = "black") +
    # scale_x_continuous(limits = c(0, 10), breaks = c(0,5,10))+
    # scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
    ggtitle(paste("DGE, ",sample,
                  "\nDN=",table(treat_gene_expression5$regu.adj)[1],
                  ", NC=",table(treat_gene_expression5$regu.adj)[2],
                  ", UP=",table(treat_gene_expression5$regu.adj)[3],
                  "\nUP/DN=",round(table(treat_gene_expression5$regu.adj)[3]/table(treat_gene_expression5$regu.adj)[1],2),
                  ", DN/UP=",round(table(treat_gene_expression5$regu.adj)[1]/table(treat_gene_expression5$regu.adj)[3],2),sep = ""))
  
}



library(cowplot)
plot_grid(plotlist=c(GEplot_MA),ncol = 4)
