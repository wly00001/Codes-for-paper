
# 700*700 .SVG

library(DESeq2)
library(ggplot2)
library(reshape)
library(tidyverse)
library(ggrepel)


gene_full_name <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/18077-13/gene_full_names + additional info.txt", header = T, sep = "\t", quote = "", na.strings = c("","NA"))
# gene_full_name <- gene_full_name[-c(1,2)]
gene_full_name <- gene_full_name %>% separate(gene_name, c("gene_name",NA),sep=",,")
gene_full_name <- gene_full_name[!duplicated(gene_full_name$gene_symbol), ]

### generating expression table from raw featurecounts outputs
id_convert <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/id convert.tbl", header = F, sep = "\t")
names(id_convert) <- c("gene_id","gene_symbol")
id_convert2 <- id_convert[!duplicated(id_convert), ]
siLUC_rep1 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep1.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep1) <- c("gene_id","siLUC_rep1")
siLUC_rep2 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep2.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep2) <- c("gene_id","siLUC_rep2")
siLUC_rep3 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep3.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep3) <- c("gene_id","siLUC_rep3")
siLUC_rep4 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep4.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep4) <- c("gene_id","siLUC_rep4")

siPCF11_rep1 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep1.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep1) <- c("gene_id","siPCF11_rep1")
siPCF11_rep2 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep2.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep2) <- c("gene_id","siPCF11_rep2")
siPCF11_rep3 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep3.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep3) <- c("gene_id","siPCF11_rep3")
siPCF11_rep4 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep4.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep4) <- c("gene_id","siPCF11_rep4")

gene_expression <- Reduce(function(x, y) merge(x, y), list(id_convert2,
                                                           gene_full_name,
                                                           siLUC_rep1,
                                                           siLUC_rep2,
                                                           siLUC_rep3,
                                                           siLUC_rep4,
                                                           siPCF11_rep1,
                                                           siPCF11_rep2,
                                                           siPCF11_rep3,
                                                           siPCF11_rep4))
gene_expression_nonMT <- subset(gene_expression, Chromosome != "MT")
gene_expression2 <- gene_expression_nonMT[-c(1,3:8)]
gene_expression3 <- aggregate(. ~ gene_symbol, gene_expression2, sum)
gene_expression4 <- subset(gene_expression3, siLUC_rep1>5&siLUC_rep2>5&siLUC_rep3>5&siLUC_rep4>5&
                             siPCF11_rep1>5&siPCF11_rep2>5&siPCF11_rep3>5&siPCF11_rep4>5)
gene_expression5 <- data.frame(gene_expression4[-c(1)],row.names=gene_expression4$gene_symbol)

colData <- data.frame("condition" = c(rep("ctrl",4),rep("siPCF11",4)),"type" = c(rep("NA",8)),row.names = names(gene_expression5))
dds <- DESeqDataSetFromMatrix(countData = gene_expression5, colData = colData, design = ~ condition)
dds <- DESeq(dds)
dds

res_siPCF11_ctrl <- results(dds, contrast=c("condition","siPCF11","ctrl"), alpha=0.05, independentFiltering=FALSE)
res_siPCF11_ctrl <- data.frame(res_siPCF11_ctrl)
res_siPCF11_ctrl_sig <- subset(res_siPCF11_ctrl, padj < 0.05 & abs(log2FoldChange)>log2(1.2))
res_siPCF11_ctrl2 <- res_siPCF11_ctrl
res_siPCF11_ctrl2$regu <- ifelse(res_siPCF11_ctrl2$padj < 0.05 & res_siPCF11_ctrl2$log2FoldChange>log2(1.2), "UP", ifelse(res_siPCF11_ctrl2$padj < 0.05 & res_siPCF11_ctrl2$log2FoldChange<(-log2(1.2)), "DN", "NO"))
res_siPCF11_ctrl2$gene_symbol <- row.names(res_siPCF11_ctrl2)

DGE_NP <- res_siPCF11_ctrl2[c(8,2)]
names(DGE_NP)[c(2)] <- "log2_NP"




List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData")

List.DGE_QSF2 <- lapply(List.DGE_QSF, function(x) { x <- x[c(1,6)]; x })
DGE_QSF <- Reduce(function(x, y) merge(x, y), List.DGE_QSF2)

DGE_Jb <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/CRISPRpas/_DGE/QSP.J_PCF11_b.v.J_PCF11_ctrl.BH.csv"), header = T)[c(1,6)]
DGE_Jd <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/CRISPRpas/_DGE/QSP.J_PCF11_d.v.J_PCF11_ctrl.BH.csv"), header = T)[c(1,6)]

DGE <- Reduce(function(x, y) merge(x, y), list(DGE_QSF,DGE_NP,DGE_Jb,DGE_Jd))


#####################
### input DE data ###

library(reshape)
library(tidyverse)
library(ggrepel)
library(cowplot)

mychart.Correlation <- function (R, method = c("pearson", "kendall", "spearman"), ...){

  x = PerformanceAnalytics::checkData(R, method = "matrix")
  if (missing(method)) method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "",
                        use = "pairwise.complete.obs",
                        method = cormeth,
                        cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) cex <- 1.5
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }

  pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)

}




tbl <- data.frame(DGE[-c(1)],row.names=DGE$gene_symbol)
tbl <- tbl %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
tbl <- na.omit(tbl)

# write.csv(tbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig3e.csv"),row.names = T)
# write.csv(tbl, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS1b.csv"),row.names = T)

# mychart.Correlation(tbl, histogram=FALSE, method = "pearson")




tbl2 <- tbl[c(1:4)]
x <- c("log2.HeLa_JTE_1.v.HeLa_DMSO","log2.HeLa_JTE_1.v.HeLa_DMSO","log2.HeLa_JTE_1.v.HeLa_DMSO",
       "log2.HeLa_JTE_10.v.HeLa_DMSO","log2.HeLa_JTE_10.v.HeLa_DMSO","log2.HepG2_JTE1.v.HepG2_DMSO")
y <- c("log2.HeLa_JTE_10.v.HeLa_DMSO","log2.HepG2_JTE1.v.HepG2_DMSO","log2.HepG2_JTE10.v.HepG2_DMSO",
       "log2.HepG2_JTE1.v.HepG2_DMSO","log2.HepG2_JTE10.v.HepG2_DMSO","log2.HepG2_JTE10.v.HepG2_DMSO")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 200) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
    
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)



cor.test(tbl$log2.HepG2_JTE10.v.HepG2_DMSO,tbl$log2_NP)
ggplot(data=tbl, aes(x=log2.HepG2_JTE10.v.HepG2_DMSO, y=log2_NP)) +
  geom_point() +
  geom_bin2d(bins = 200) +
  theme_bw() +
  coord_fixed(ratio = 1)+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_vline(xintercept=0, col="black", linetype="solid") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  labs(x = "HepG2_JTE10", y = "log2_NP")


cor.test(tbl$log2.HeLa_JTE_10.v.HeLa_DMSO,tbl$log2_NP)
ggplot(data=tbl, aes(x=log2.HeLa_JTE_10.v.HeLa_DMSO, y=log2_NP)) +
  geom_point() +
  geom_bin2d(bins = 250) +
  theme_bw() +
  coord_fixed(ratio = 1)+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_vline(xintercept=0, col="black", linetype="solid") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  labs(x = "HeLa_JTE10", y = "log2_NP")











# 700*700 .SVG

library(DESeq2)
library(ggplot2)
library(reshape)
library(tidyverse)
library(ggrepel)


gene_full_name <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/18077-13/gene_full_names + additional info.txt", header = T, sep = "\t", quote = "", na.strings = c("","NA"))
# gene_full_name <- gene_full_name[-c(1,2)]
gene_full_name <- gene_full_name %>% separate(gene_name, c("gene_name",NA),sep=",,")
gene_full_name <- gene_full_name[!duplicated(gene_full_name$gene_symbol), ]

### generating expression table from raw featurecounts outputs
id_convert <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/id convert.tbl", header = F, sep = "\t")
names(id_convert) <- c("gene_id","gene_symbol")
id_convert2 <- id_convert[!duplicated(id_convert), ]
siLUC_rep1 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep1.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep1) <- c("gene_id","siLUC_rep1")
siLUC_rep2 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep2.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep2) <- c("gene_id","siLUC_rep2")
siLUC_rep3 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep3.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep3) <- c("gene_id","siLUC_rep3")
siLUC_rep4 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siLUC_rep4.featurecounts.Rmatrix.txt", header = T)
names(siLUC_rep4) <- c("gene_id","siLUC_rep4")

siPCF11_rep1 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep1.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep1) <- c("gene_id","siPCF11_rep1")
siPCF11_rep2 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep2.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep2) <- c("gene_id","siPCF11_rep2")
siPCF11_rep3 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep3.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep3) <- c("gene_id","siPCF11_rep3")
siPCF11_rep4 <- read.table("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/featureCounts/3mRNAseq_siPCF11_rep4.featurecounts.Rmatrix.txt", header = T)
names(siPCF11_rep4) <- c("gene_id","siPCF11_rep4")

gene_expression <- Reduce(function(x, y) merge(x, y), list(id_convert2,
                                                           gene_full_name,
                                                           siLUC_rep1,
                                                           siLUC_rep2,
                                                           siLUC_rep3,
                                                           siLUC_rep4,
                                                           siPCF11_rep1,
                                                           siPCF11_rep2,
                                                           siPCF11_rep3,
                                                           siPCF11_rep4))
gene_expression_nonMT <- subset(gene_expression, Chromosome != "MT")
gene_expression2 <- gene_expression_nonMT[-c(1,3:8)]
gene_expression3 <- aggregate(. ~ gene_symbol, gene_expression2, sum)
gene_expression4 <- subset(gene_expression3, siLUC_rep1>5&siLUC_rep2>5&siLUC_rep3>5&siLUC_rep4>5&
                             siPCF11_rep1>5&siPCF11_rep2>5&siPCF11_rep3>5&siPCF11_rep4>5)
gene_expression5 <- data.frame(gene_expression4[-c(1)],row.names=gene_expression4$gene_symbol)

colData <- data.frame("condition" = c(rep("ctrl",4),rep("siPCF11",4)),"type" = c(rep("NA",8)),row.names = names(gene_expression5))
dds <- DESeqDataSetFromMatrix(countData = gene_expression5, colData = colData, design = ~ condition)
dds <- DESeq(dds)
dds

res_siPCF11_ctrl <- results(dds, contrast=c("condition","siPCF11","ctrl"), alpha=0.05, independentFiltering=FALSE)
res_siPCF11_ctrl <- data.frame(res_siPCF11_ctrl)
res_siPCF11_ctrl_sig <- subset(res_siPCF11_ctrl, padj < 0.05 & abs(log2FoldChange)>log2(1.2))
res_siPCF11_ctrl2 <- res_siPCF11_ctrl
res_siPCF11_ctrl2$regu <- ifelse(res_siPCF11_ctrl2$padj < 0.05 & res_siPCF11_ctrl2$log2FoldChange>log2(1.2), "UP", ifelse(res_siPCF11_ctrl2$padj < 0.05 & res_siPCF11_ctrl2$log2FoldChange<(-log2(1.2)), "DN", "NO"))
res_siPCF11_ctrl2$gene_symbol <- row.names(res_siPCF11_ctrl2)

DGE_NP <- res_siPCF11_ctrl2[c(8,2)]
names(DGE_NP)[c(2)] <- "log2_NP"




List.DGE_QSF <- readRDS("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-08/_DGE_ChrRNA/List.DGE_QSF_1_10.RData")

filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-09/_DGE/FC1.2/", pattern="*.csv", full.names=TRUE)
filenames
ldf <- lapply(filenames, read.csv)
names(ldf) <- substr(filenames,66,nchar(filenames)-7)
names(ldf) <- gsub("YC_103_pA_1_1","24DnJn",names(ldf))
names(ldf) <- gsub("YC_103_pA_1_2","24DnJp",names(ldf))
names(ldf) <- gsub("YC_103_pA_1_3","24DpJn",names(ldf))
names(ldf) <- gsub("YC_103_pA_1_4","24DpJp",names(ldf))
names(ldf) <- gsub("YC_103_pA_2_5","48DnJn",names(ldf))
names(ldf) <- gsub("YC_103_pA_2_6","48DnJp",names(ldf))
names(ldf) <- gsub("YC_103_pA_2_7","48DpJn",names(ldf))
names(ldf) <- gsub("YC_103_pA_2_8","48DpJp",names(ldf))
# ldf <- lapply(ldf, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
ldf <- lapply(ldf, function(x) { x[sapply(x, is.infinite)] <- NA; x })
ldf <- lapply(ldf, function(x) { x<-na.omit(x); x })

new.ldf <- list()
for (n in 1:length(ldf)) {
  
  name <- names(ldf)[c(n)]
  
  df <- ldf[[n]]
  
  names(df) <- gsub("YC_103_pA_1_1","24DnJn",names(df))
  names(df) <- gsub("YC_103_pA_1_2","24DnJp",names(df))
  names(df) <- gsub("YC_103_pA_1_3","24DpJn",names(df))
  names(df) <- gsub("YC_103_pA_1_4","24DpJp",names(df))
  names(df) <- gsub("YC_103_pA_2_5","48DnJn",names(df))
  names(df) <- gsub("YC_103_pA_2_6","48DnJp",names(df))
  names(df) <- gsub("YC_103_pA_2_7","48DpJn",names(df))
  names(df) <- gsub("YC_103_pA_2_8","48DpJp",names(df))
  
  new.ldf[[name]] <- df
  
  
  
}

filenames <- list.files("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/siRNA1204/_DGE/", pattern="*.csv", full.names=TRUE)
filenames
ldf <- lapply(filenames, read.csv)
names(ldf) <- substr(filenames,69,nchar(filenames)-7)
# ldf <- lapply(ldf, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
ldf <- lapply(ldf, function(x) { x[sapply(x, is.infinite)] <- NA; x })
ldf <- lapply(ldf, function(x) { x<-na.omit(x); x })
new.ldf2 <- ldf

List.DGE_QSF2 <- lapply(c(List.DGE_QSF,new.ldf2,new.ldf[c(1,6,2,7,5,10)]), function(x) { x <- x[c(1,6)]; x })
DGE_QSF <- Reduce(function(x, y) merge(x, y), List.DGE_QSF2[-c(7,9,10)])

# DGE_Jb <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/CRISPRpas/_DGE/QSP.J_PCF11_b.v.J_PCF11_ctrl.BH.csv"), header = T)[c(1,6)]
# DGE_Jd <- read.csv(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/CRISPRpas/_DGE/QSP.J_PCF11_d.v.J_PCF11_ctrl.BH.csv"), header = T)[c(1,6)]
# DGE <- Reduce(function(x, y) merge(x, y), list(DGE_QSF,DGE_NP,DGE_Jb,DGE_Jd))

DGE <- Reduce(function(x, y) merge(x, y), list(DGE_QSF,DGE_NP))





tbl <- data.frame(DGE[-c(1)],row.names=DGE$gene_symbol)
tbl <- tbl %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
tbl <- na.omit(tbl)



tbl2 <- tbl[c(3:6)]
# write.csv(tbl2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS6b.csv"),row.names = T)

x <- c("log2.HepG2_JTE1.v.HepG2_DMSO","log2.HepG2_JTE1.v.HepG2_DMSO","log2.HepG2_JTE1.v.HepG2_DMSO",
       "log2.HepG2_JTE10.v.HepG2_DMSO","log2.HepG2_JTE10.v.HepG2_DMSO",
       "log2.U937_JTE607_6h.v.U937_Mock_6h")
y <- c("log2.HepG2_JTE10.v.HepG2_DMSO","log2.U937_JTE607_6h.v.U937_Mock_6h","log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "log2.U937_JTE607_6h.v.U937_Mock_6h","log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 200) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-3, 3))+
    scale_y_continuous(limits = c(-3, 3))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)















