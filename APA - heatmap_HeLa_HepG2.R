workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"

#######  APA heatmap

library(reshape2)
library(tidyverse)
library(RColorBrewer)

#REDu

###input generated tbl
list.REDu <- readRDS(paste0(workplace,"/linux/lwang/allJTE_06_07_sum0909/result/HeLa_HepG2/List.REDu.RData"))[c(1:4)]
names(list.REDu) <- c("HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h")

#REDi

###input generated tbl
list.REDi <- readRDS(paste0(workplace,"/linux/lwang/allJTE_06_07_sum0909/result/HeLa_HepG2/List.REDi_1vAll.RData"))[c(1:4)]
names(list.REDi) <- c("HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h")

list.REDu <- lapply(list.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDu <- lapply(list.REDu, function(x) { x<-na.omit(x); x })
list.REDi <- lapply(list.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDi <- lapply(list.REDi, function(x) { x<-na.omit(x); x })

list.REDu2 <- list()
list.REDi2 <- list()

for (sample in names(list.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- list.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  list.REDu2[[sample]] <- df
  
}

for (sample in names(list.REDi)) {
  
  # sample="HeLa_JTE_01"
  
  df <- list.REDi[[sample]]
  names(df)[c(4:13)] <- paste0(names(df)[c(4:13)],".",sample)
  list.REDi2[[sample]] <- df
  
}

List.REDu.8h3 <- lapply(list.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3)
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO" | heatmapTbl.REDu[[7]]!="NO" | heatmapTbl.REDu[[9]]!="NO")
heatmapTbl.REDu2.HeLa <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO")
heatmapTbl.REDu2.HepG2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[7]]!="NO" | heatmapTbl.REDu[[9]]!="NO")

heatmapTbl.REDu3 <- data.frame(heatmapTbl.REDu2[-c(1,3,5,7,9)],row.names=heatmapTbl.REDu2$gene_symbol)
# heatmapTbl.REDu3 <- round(heatmapTbl.REDu3, 2)
heatmapTbl.REDu3 <- heatmapTbl.REDu3[Reduce(`&`, lapply(heatmapTbl.REDu3, is.finite)),]
names(heatmapTbl.REDu3) <- gsub("REDu1.","",names(heatmapTbl.REDu3))


heatmapTbl.REDu4 <- as.matrix(heatmapTbl.REDu3)

# # heatmap with clustering using the euclidean distance as dissimilarity measure
# heatmap.2(df3, col=bluered(100),
#           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
#           scale="none", breaks = seq(-4, 4, length.out = 101))

# heatmap with clustering using Pearson correlation as dissimilarity measures

# Pairwise correlation between samples (columns)
cols.cor <- cor(heatmapTbl.REDu4, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(heatmapTbl.REDu4), use = "pairwise.complete.obs", method = "pearson")

## Row- and column-wise clustering using correlation
hclust.col <- hclust(as.dist(1-cols.cor))
hclust.row <- hclust(as.dist(1-rows.cor))

# Plot the heatmap
library(gplots)
hm.REDu <- heatmap.2(heatmapTbl.REDu4, col=bluered(100),
                     density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
                     scale="none", breaks = seq(-5, 5, length.out = 101),
                     # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
                     Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("all, N = ",nrow(heatmapTbl.REDu4)))
dev.off()

# write.csv(heatmapTbl.REDu2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5b.csv"),row.names = T)




# library(cmstatr)
# heatmapTbl.REDu3$CV <- apply(heatmapTbl.REDu3, 1, cv)
# 
# heatmapTbl.REDu3.highVar <- subset(heatmapTbl.REDu3, abs(CV) > 1)[-5]
# heatmapTbl.REDu3.lowVar <- subset(heatmapTbl.REDu3, abs(CV) <= 1)[-5]
# 
# 
# 
# heatmapTbl.REDu4 <- as.matrix(heatmapTbl.REDu3.highVar)
# 
# # # heatmap with clustering using the euclidean distance as dissimilarity measure
# # heatmap.2(df3, col=bluered(100),
# #           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
# #           scale="none", breaks = seq(-4, 4, length.out = 101))
# 
# # heatmap with clustering using Pearson correlation as dissimilarity measures
# 
# # Pairwise correlation between samples (columns)
# cols.cor <- cor(heatmapTbl.REDu4, use = "pairwise.complete.obs", method = "pearson")
# # Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatmapTbl.REDu4), use = "pairwise.complete.obs", method = "pearson")
# 
# ## Row- and column-wise clustering using correlation
# hclust.col <- hclust(as.dist(1-cols.cor))
# hclust.row <- hclust(as.dist(1-rows.cor))
# 
# # Plot the heatmap
# library(gplots)
# hm.REDu <- heatmap.2(heatmapTbl.REDu4, col=bluered(100),
#                      density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
#                      scale="none", breaks = seq(-5, 5, length.out = 101),
#                      # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
#                      Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("CV > 1, N = ",nrow(heatmapTbl.REDu4)))
# dev.off()
# 
# 
# 
# heatmapTbl.REDu4 <- as.matrix(heatmapTbl.REDu3.lowVar)
# 
# # # heatmap with clustering using the euclidean distance as dissimilarity measure
# # heatmap.2(df3, col=bluered(100),
# #           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
# #           scale="none", breaks = seq(-4, 4, length.out = 101))
# 
# # heatmap with clustering using Pearson correlation as dissimilarity measures
# 
# # Pairwise correlation between samples (columns)
# cols.cor <- cor(heatmapTbl.REDu4, use = "pairwise.complete.obs", method = "pearson")
# # Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatmapTbl.REDu4), use = "pairwise.complete.obs", method = "pearson")
# 
# ## Row- and column-wise clustering using correlation
# hclust.col <- hclust(as.dist(1-cols.cor))
# hclust.row <- hclust(as.dist(1-rows.cor))
# 
# # Plot the heatmap
# library(gplots)
# hm.REDu <- heatmap.2(heatmapTbl.REDu4, col=bluered(100),
#                      density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
#                      scale="none", breaks = seq(-5, 5, length.out = 101),
#                      # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
#                      Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("CV <= 1, N = ",nrow(heatmapTbl.REDu4)))
# dev.off()



## gene selection, select genes that are consistently regulated in each cell line: a) two concentration data show the same direction; b) 10uM shows a larger change than 1uM (either direction). 


heatmapDF.REDu <- heatmapTbl.REDu2[-c(3,5,7,9)]
names(heatmapDF.REDu) <- c("gene_symbol","HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h")

heatmapDF.REDu_HeLa <- heatmapDF.REDu[c(1,2,3)]
heatmapDF.REDu_HepG2 <- heatmapDF.REDu[c(1,4,5)]

heatmapDF.REDu_HeLa_le <- subset(heatmapDF.REDu_HeLa, heatmapDF.REDu_HeLa[[2]]>0 & heatmapDF.REDu_HeLa[[3]]>0 & heatmapDF.REDu_HeLa[[2]]<heatmapDF.REDu_HeLa[[3]])
heatmapDF.REDu_HeLa_sh <- subset(heatmapDF.REDu_HeLa, heatmapDF.REDu_HeLa[[2]]<0 & heatmapDF.REDu_HeLa[[3]]<0 & heatmapDF.REDu_HeLa[[2]]>heatmapDF.REDu_HeLa[[3]])

heatmapDF.REDu_HepG2_le <- subset(heatmapDF.REDu_HepG2, heatmapDF.REDu_HepG2[[2]]>0 & heatmapDF.REDu_HepG2[[3]]>0 & heatmapDF.REDu_HepG2[[2]]<heatmapDF.REDu_HepG2[[3]])
heatmapDF.REDu_HepG2_sh <- subset(heatmapDF.REDu_HepG2, heatmapDF.REDu_HepG2[[2]]<0 & heatmapDF.REDu_HepG2[[3]]<0 & heatmapDF.REDu_HepG2[[2]]>heatmapDF.REDu_HepG2[[3]])

# library(UpSetR)
# LElist <- list(LE_HeLa = heatmapDF.REDu_HeLa_le$gene_symbol,
#                LE_HepG2 = heatmapDF.REDu_HepG2_le$gene_symbol)
# upset(fromList(LElist), nsets = 2, nintersects = NA, sets = names(LElist), keep.order = T, order.by = "freq", text.scale = 1.5)
# 
# SHlist <- list(SH_HeLa = heatmapDF.REDu_HeLa_sh$gene_symbol,
#                SH_HepG2 = heatmapDF.REDu_HepG2_sh$gene_symbol)
# upset(fromList(SHlist), nsets = 2, nintersects = NA, sets = names(SHlist), keep.order = T, order.by = "freq", text.scale = 1.5)


setList <- list()

setList[["all"]] <- data.frame("gene_symbol"=heatmapTbl.REDu$gene_symbol)

setList[["LE_common"]] <- data.frame("gene_symbol"=intersect(heatmapDF.REDu_HeLa_le$gene_symbol,heatmapDF.REDu_HepG2_le$gene_symbol))
LE_HeLa <- data.frame("gene_symbol"=setdiff(heatmapDF.REDu_HeLa_le$gene_symbol,heatmapDF.REDu_HepG2_le$gene_symbol))
LE_HeLa <- subset(LE_HeLa, gene_symbol %in% heatmapTbl.REDu2.HeLa$gene_symbol)
setList[["LE_HeLa"]] <- LE_HeLa
LE_HepG2 <- data.frame("gene_symbol"=setdiff(heatmapDF.REDu_HepG2_le$gene_symbol,heatmapDF.REDu_HeLa_le$gene_symbol))
LE_HepG2 <- subset(LE_HepG2, gene_symbol %in% heatmapTbl.REDu2.HepG2$gene_symbol)
setList[["LE_HepG2"]] <- LE_HepG2

setList[["SH_common"]] <- data.frame("gene_symbol"=intersect(heatmapDF.REDu_HeLa_sh$gene_symbol,heatmapDF.REDu_HepG2_sh$gene_symbol))
SH_HeLa <- data.frame("gene_symbol"=setdiff(heatmapDF.REDu_HeLa_sh$gene_symbol,heatmapDF.REDu_HepG2_sh$gene_symbol))
SH_HeLa <- subset(SH_HeLa, gene_symbol %in% heatmapTbl.REDu2.HeLa$gene_symbol)
setList[["SH_HeLa"]] <- SH_HeLa
SH_HepG2 <- data.frame("gene_symbol"=setdiff(heatmapDF.REDu_HepG2_sh$gene_symbol,heatmapDF.REDu_HeLa_sh$gene_symbol))
SH_HepG2 <- subset(SH_HepG2, gene_symbol %in% heatmapTbl.REDu2.HepG2$gene_symbol)
setList[["SH_HepG2"]] <- SH_HepG2


library(UpSetR)
LElist <- list(LE_HeLa = c(setList[["LE_common"]]$gene_symbol,setList[["LE_HeLa"]]$gene_symbol),
               LE_HepG2 = c(setList[["LE_common"]]$gene_symbol,setList[["LE_HepG2"]]$gene_symbol))
upset(fromList(LElist), nsets = 2, nintersects = NA, sets = names(LElist), keep.order = T, order.by = "freq", text.scale = 1.5)

SHlist <- list(SH_HeLa = c(setList[["SH_common"]]$gene_symbol,setList[["SH_HeLa"]]$gene_symbol),
               SH_HepG2 = c(setList[["SH_common"]]$gene_symbol,setList[["SH_HepG2"]]$gene_symbol))
upset(fromList(SHlist), nsets = 2, nintersects = NA, sets = names(SHlist), keep.order = T, order.by = "freq", text.scale = 1.5)



# saveRDS(setList, file=paste0("./setList_REDu_HeLa_HepG2.RData"))




library(UpSetR)
list.LE_HeLa.SH_HepG2 <- list(LE_HeLa = heatmapDF.REDu_HeLa_le$gene_symbol,
                              SH_HepG2 = heatmapDF.REDu_HepG2_sh$gene_symbol)
upset(fromList(list.LE_HeLa.SH_HepG2), nsets = 2, nintersects = NA, sets = names(list.LE_HeLa.SH_HepG2), keep.order = T, order.by = "freq", text.scale = 1.5)

list.SH_HeLa.LE_HepG2 <- list(SH_HeLa = heatmapDF.REDu_HeLa_sh$gene_symbol,
                              LE_HepG2 = heatmapDF.REDu_HepG2_le$gene_symbol)
upset(fromList(list.SH_HeLa.LE_HepG2), nsets = 2, nintersects = NA, sets = names(list.SH_HeLa.LE_HepG2), keep.order = T, order.by = "freq", text.scale = 1.5)















List.REDi.8h3 <- lapply(list.REDi2, function(x) { x <- x[c(1,8,13)]; x })
heatmapTbl.REDi <- Reduce(function(x, y) merge(x, y), List.REDi.8h3)
heatmapTbl.REDi2 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO" | heatmapTbl.REDi[[5]]!="NO" | heatmapTbl.REDi[[7]]!="NO" | heatmapTbl.REDi[[9]]!="NO")
heatmapTbl.REDi2.HeLa <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO" | heatmapTbl.REDi[[5]]!="NO")
heatmapTbl.REDi2.HepG2 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[7]]!="NO" | heatmapTbl.REDi[[9]]!="NO")

heatmapTbl.REDi3 <- data.frame(heatmapTbl.REDi2[-c(1,3,5,7,9)],row.names=heatmapTbl.REDi2$gene_symbol)
# heatmapTbl.REDi3 <- round(heatmapTbl.REDi3, 2)
heatmapTbl.REDi3 <- heatmapTbl.REDi3[Reduce(`&`, lapply(heatmapTbl.REDi3, is.finite)),]
names(heatmapTbl.REDi3) <- gsub("REDi1.","",names(heatmapTbl.REDi3))


heatmapTbl.REDi4 <- as.matrix(heatmapTbl.REDi3)

# # heatmap with clustering using the euclidean distance as dissimilarity measure
# heatmap.2(df3, col=bluered(100),
#           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
#           scale="none", breaks = seq(-4, 4, length.out = 101))

# heatmap with clustering using Pearson correlation as dissimilarity measures

# Pairwise correlation between samples (columns)
cols.cor <- cor(heatmapTbl.REDi4, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(heatmapTbl.REDi4), use = "pairwise.complete.obs", method = "pearson")

## Row- and column-wise clustering using correlation
hclust.col <- hclust(as.dist(1-cols.cor))
hclust.row <- hclust(as.dist(1-rows.cor))

# Plot the heatmap
library(gplots)
hm.REDi <- heatmap.2(heatmapTbl.REDi4, col=bluered(100),
                     density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
                     scale="none", breaks = seq(-5, 5, length.out = 101),
                     # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
                     Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("all, N = ",nrow(heatmapTbl.REDi4)))
dev.off()

# write.csv(heatmapTbl.REDi2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig5gh.csv"),row.names = T)




# library(cmstatr)
# heatmapTbl.REDi3$CV <- apply(heatmapTbl.REDi3, 1, cv)
# 
# heatmapTbl.REDi3.highVar <- subset(heatmapTbl.REDi3, abs(CV) > 1)[-5]
# heatmapTbl.REDi3.lowVar <- subset(heatmapTbl.REDi3, abs(CV) <= 1)[-5]
# 
# 
# 
# heatmapTbl.REDi4 <- as.matrix(heatmapTbl.REDi3.highVar)
# 
# # # heatmap with clustering using the euclidean distance as dissimilarity measure
# # heatmap.2(df3, col=bluered(100),
# #           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
# #           scale="none", breaks = seq(-4, 4, length.out = 101))
# 
# # heatmap with clustering using Pearson correlation as dissimilarity measures
# 
# # Pairwise correlation between samples (columns)
# cols.cor <- cor(heatmapTbl.REDi4, use = "pairwise.complete.obs", method = "pearson")
# # Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatmapTbl.REDi4), use = "pairwise.complete.obs", method = "pearson")
# 
# ## Row- and column-wise clustering using correlation
# hclust.col <- hclust(as.dist(1-cols.cor))
# hclust.row <- hclust(as.dist(1-rows.cor))
# 
# # Plot the heatmap
# library(gplots)
# hm.REDi <- heatmap.2(heatmapTbl.REDi4, col=bluered(100),
#                      density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
#                      scale="none", breaks = seq(-5, 5, length.out = 101),
#                      # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
#                      Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("CV > 1, N = ",nrow(heatmapTbl.REDi4)))
# dev.off()
# 
# 
# 
# heatmapTbl.REDi4 <- as.matrix(heatmapTbl.REDi3.lowVar)
# 
# # # heatmap with clustering using the euclidean distance as dissimilarity measure
# # heatmap.2(df3, col=bluered(100),
# #           density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 1,
# #           scale="none", breaks = seq(-4, 4, length.out = 101))
# 
# # heatmap with clustering using Pearson correlation as dissimilarity measures
# 
# # Pairwise correlation between samples (columns)
# cols.cor <- cor(heatmapTbl.REDi4, use = "pairwise.complete.obs", method = "pearson")
# # Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatmapTbl.REDi4), use = "pairwise.complete.obs", method = "pearson")
# 
# ## Row- and column-wise clustering using correlation
# hclust.col <- hclust(as.dist(1-cols.cor))
# hclust.row <- hclust(as.dist(1-rows.cor))
# 
# # Plot the heatmap
# library(gplots)
# hm.REDi <- heatmap.2(heatmapTbl.REDi4, col=bluered(100),
#                      density.info = "none", trace="none", na.color="gray", lhei = c(2, 10), labRow = F, cexCol = 0.6,
#                      scale="none", breaks = seq(-5, 5, length.out = 101),
#                      # Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row))
#                      Colv = as.dendrogram(hclust.col),Rowv = as.dendrogram(hclust.row),main = paste0("CV <= 1, N = ",nrow(heatmapTbl.REDi4)))
# dev.off()


## gene selection, select genes that are consistently regulated in each cell line: a) two concentration data show the same direction; b) 10uM shows a larger change than 1uM (either direction). 


heatmapDF.REDi <- heatmapTbl.REDi2[-c(3,5,7,9)]
names(heatmapDF.REDi) <- c("gene_symbol","HeLa_1_9h","HeLa_10_9h","HepG2_1_7h","HepG2_10_7h")

heatmapDF.REDi_HeLa <- heatmapDF.REDi[c(1,2,3)]
heatmapDF.REDi_HepG2 <- heatmapDF.REDi[c(1,4,5)]

heatmapDF.REDi_HeLa_su <- subset(heatmapDF.REDi_HeLa, heatmapDF.REDi_HeLa[[2]]>0 & heatmapDF.REDi_HeLa[[3]]>0 & heatmapDF.REDi_HeLa[[2]]<heatmapDF.REDi_HeLa[[3]])
heatmapDF.REDi_HeLa_ac <- subset(heatmapDF.REDi_HeLa, heatmapDF.REDi_HeLa[[2]]<0 & heatmapDF.REDi_HeLa[[3]]<0 & heatmapDF.REDi_HeLa[[2]]>heatmapDF.REDi_HeLa[[3]])

heatmapDF.REDi_HepG2_su <- subset(heatmapDF.REDi_HepG2, heatmapDF.REDi_HepG2[[2]]>0 & heatmapDF.REDi_HepG2[[3]]>0 & heatmapDF.REDi_HepG2[[2]]<heatmapDF.REDi_HepG2[[3]])
heatmapDF.REDi_HepG2_ac <- subset(heatmapDF.REDi_HepG2, heatmapDF.REDi_HepG2[[2]]<0 & heatmapDF.REDi_HepG2[[3]]<0 & heatmapDF.REDi_HepG2[[2]]>heatmapDF.REDi_HepG2[[3]])

# library(UpSetR)
# SUlist <- list(SU_HeLa = heatmapDF.REDi_HeLa_su$gene_symbol,
#                SU_HepG2 = heatmapDF.REDi_HepG2_su$gene_symbol)
# upset(fromList(SUlist), nsets = 2, nintersects = NA, sets = names(SUlist), keep.order = T, order.by = "freq", text.scale = 1.5)
# 
# AClist <- list(AC_HeLa = heatmapDF.REDi_HeLa_ac$gene_symbol,
#                AC_HepG2 = heatmapDF.REDi_HepG2_ac$gene_symbol)
# upset(fromList(AClist), nsets = 2, nintersects = NA, sets = names(AClist), keep.order = T, order.by = "freq", text.scale = 1.5)


setList <- list()

setList[["all"]] <- data.frame("gene_symbol"=heatmapTbl.REDi$gene_symbol)

setList[["SU_common"]] <- data.frame("gene_symbol"=intersect(heatmapDF.REDi_HeLa_su$gene_symbol,heatmapDF.REDi_HepG2_su$gene_symbol))
SU_HeLa <- data.frame("gene_symbol"=setdiff(heatmapDF.REDi_HeLa_su$gene_symbol,heatmapDF.REDi_HepG2_su$gene_symbol))
SU_HeLa <- subset(SU_HeLa, gene_symbol %in% heatmapTbl.REDi2.HeLa$gene_symbol)
setList[["SU_HeLa"]] <- SU_HeLa
SU_HepG2 <- data.frame("gene_symbol"=setdiff(heatmapDF.REDi_HepG2_su$gene_symbol,heatmapDF.REDi_HeLa_su$gene_symbol))
SU_HepG2 <- subset(SU_HepG2, gene_symbol %in% heatmapTbl.REDi2.HepG2$gene_symbol)
setList[["SU_HepG2"]] <- SU_HepG2

setList[["AC_common"]] <- data.frame("gene_symbol"=intersect(heatmapDF.REDi_HeLa_ac$gene_symbol,heatmapDF.REDi_HepG2_ac$gene_symbol))
AC_HeLa <- data.frame("gene_symbol"=setdiff(heatmapDF.REDi_HeLa_ac$gene_symbol,heatmapDF.REDi_HepG2_ac$gene_symbol))
AC_HeLa <- subset(AC_HeLa, gene_symbol %in% heatmapTbl.REDi2.HeLa$gene_symbol)
setList[["AC_HeLa"]] <- AC_HeLa
AC_HepG2 <- data.frame("gene_symbol"=setdiff(heatmapDF.REDi_HepG2_ac$gene_symbol,heatmapDF.REDi_HeLa_ac$gene_symbol))
AC_HepG2 <- subset(AC_HepG2, gene_symbol %in% heatmapTbl.REDi2.HepG2$gene_symbol)
setList[["AC_HepG2"]] <- AC_HepG2


library(UpSetR)
SUlist <- list(SU_HeLa = c(setList[["SU_common"]]$gene_symbol,setList[["SU_HeLa"]]$gene_symbol),
               SU_HepG2 = c(setList[["SU_common"]]$gene_symbol,setList[["SU_HepG2"]]$gene_symbol))
upset(fromList(SUlist), nsets = 2, nintersects = NA, sets = names(SUlist), keep.order = T, order.by = "freq", text.scale = 1.5)

AClist <- list(AC_HeLa = c(setList[["AC_common"]]$gene_symbol,setList[["AC_HeLa"]]$gene_symbol),
               AC_HepG2 = c(setList[["AC_common"]]$gene_symbol,setList[["AC_HepG2"]]$gene_symbol))
upset(fromList(AClist), nsets = 2, nintersects = NA, sets = names(AClist), keep.order = T, order.by = "freq", text.scale = 1.5)



# saveRDS(setList, file=paste0("./setList_REDi_HeLa_HepG2.RData"))






library(UpSetR)
list.SU_HeLa.AC_HepG2 <- list(SU_HeLa = heatmapDF.REDi_HeLa_su$gene_symbol,
                              AC_HepG2 = heatmapDF.REDi_HepG2_ac$gene_symbol)
upset(fromList(list.SU_HeLa.AC_HepG2), nsets = 2, nintersects = NA, sets = names(list.SU_HeLa.AC_HepG2), keep.order = T, order.by = "freq", text.scale = 1.5)

list.AC_HeLa.SU_HepG2 <- list(AC_HeLa = heatmapDF.REDi_HeLa_ac$gene_symbol,
                              SU_HepG2 = heatmapDF.REDi_HepG2_su$gene_symbol)
upset(fromList(list.AC_HeLa.SU_HepG2), nsets = 2, nintersects = NA, sets = names(list.AC_HeLa.SU_HepG2), keep.order = T, order.by = "freq", text.scale = 1.5)





