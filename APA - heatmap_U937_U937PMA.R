workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"


#######  APA heatmap, 7h

library(reshape2)
library(tidyverse)
library(RColorBrewer)

#REDu

###input generated tbl
list.REDu <- readRDS(paste0(workplace,"/linux/lwang/JTEonly_APA_strPAS/U937_1h6h/result/U937_1h6h/List.REDu.RData"))[c(4,5)]
# list.REDu <- readRDS(paste0(workplace,"/linux/lwang/allJTE_06_07_sum0313/U937_read2/result/U937_1h6h/List.REDu.RData"))[c(4,5)]
names(list.REDu) <- c("U937_1_7h","U937PMA_1_7h")

#REDi

###input generated tbl
list.REDi <- readRDS(paste0(workplace,"/linux/lwang/JTEonly_APA_strPAS/U937_1h6h/result/U937_1h6h/List.REDi_1vAll.RData"))[c(4,5)]
# list.REDi <- readRDS(paste0(workplace,"/linux/lwang/allJTE_06_07_sum0313/U937_read2/result/U937_1h6h/List.REDi_1vAll.RData"))[c(4,5)]
names(list.REDi) <- c("U937_1_7h","U937PMA_1_7h")

list.REDu <- lapply(list.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDu <- lapply(list.REDu, function(x) { x<-na.omit(x); x })
list.REDi <- lapply(list.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDi <- lapply(list.REDi, function(x) { x<-na.omit(x); x })

list.REDu2 <- list()
list.REDi2 <- list()

for (sample in names(list.REDu)) {
  
  # sample="U937_JTE_01"
  
  df <- list.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  list.REDu2[[sample]] <- df
  
}

for (sample in names(list.REDi)) {
  
  # sample="U937_JTE_01"
  
  df <- list.REDi[[sample]]
  names(df)[c(4:13)] <- paste0(names(df)[c(4:13)],".",sample)
  list.REDi2[[sample]] <- df
  
}

List.REDu.8h3 <- lapply(list.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3)
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO")

heatmapTbl.REDu2.U937 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO")
heatmapTbl.REDu2.U937PMA <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[5]]!="NO")

library(UpSetR)
list.U937_U937PMA <- list(sig.LE_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="lengthened")$gene_symbol,
                          sig.LE_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="lengthened")$gene_symbol,
                          nonsig.LE_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]>0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.LE_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]>0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol,
                          sig.SH_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="shortened")$gene_symbol,
                          sig.SH_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="shortened")$gene_symbol,
                          nonsig.SH_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]<0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.SH_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]<0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol)
upset(fromList(list.U937_U937PMA), nsets = 8, nintersects = NA, sets = names(list.U937_U937PMA), keep.order = T, order.by = "freq", text.scale = 1.5)

# write.csv(heatmapTbl.REDu2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig7g.csv"),row.names = T)




List.REDi.8h3 <- lapply(list.REDi2, function(x) { x <- x[c(1,8,13)]; x })
heatmapTbl.REDi <- Reduce(function(x, y) merge(x, y), List.REDi.8h3)
heatmapTbl.REDi2 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO" | heatmapTbl.REDi[[5]]!="NO")

heatmapTbl.REDi2.U937 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO")
heatmapTbl.REDi2.U937PMA <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[5]]!="NO")

library(UpSetR)
list.U937_U937PMA <- list(sig.SU_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="lengthened")$gene_symbol,
                          sig.SU_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="lengthened")$gene_symbol,
                          nonsig.SU_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]>0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.SU_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]>0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol,
                          sig.AC_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="shortened")$gene_symbol,
                          sig.AC_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="shortened")$gene_symbol,
                          nonsig.AC_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]<0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.AC_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]<0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol)
upset(fromList(list.U937_U937PMA), nsets = 8, nintersects = NA, sets = names(list.U937_U937PMA), keep.order = T, order.by = "freq", text.scale = 1.5)









#######  APA heatmap, 2h

library(reshape2)
library(tidyverse)
library(RColorBrewer)

#REDu

###input generated tbl
list.REDu <- readRDS(paste0(workplace,"/linux/lwang/JTEonly_APA_strPAS/U937_1h6h/result/U937_1h6h/List.REDu.RData"))[c(1,2)]
names(list.REDu) <- c("U937_1_2h","U937PMA_1_2h")

#REDi

###input generated tbl
list.REDi <- readRDS(paste0(workplace,"/linux/lwang/JTEonly_APA_strPAS/U937_1h6h/result/U937_1h6h/List.REDi.RData"))[c(1,2)]
names(list.REDi) <- c("U937_1_2h","U937PMA_1_2h")

list.REDu <- lapply(list.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDu <- lapply(list.REDu, function(x) { x<-na.omit(x); x })
list.REDi <- lapply(list.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
list.REDi <- lapply(list.REDi, function(x) { x<-na.omit(x); x })

list.REDu2 <- list()
list.REDi2 <- list()

for (sample in names(list.REDu)) {
  
  # sample="U937_JTE_01"
  
  df <- list.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  list.REDu2[[sample]] <- df
  
}

for (sample in names(list.REDi)) {
  
  # sample="U937_JTE_01"
  
  df <- list.REDi[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  list.REDi2[[sample]] <- df
  
}

List.REDu.8h3 <- lapply(list.REDu2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDu <- Reduce(function(x, y) merge(x, y), List.REDu.8h3)
heatmapTbl.REDu2 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO" | heatmapTbl.REDu[[5]]!="NO")

heatmapTbl.REDu2.U937 <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[3]]!="NO")
heatmapTbl.REDu2.U937PMA <- subset(heatmapTbl.REDu, heatmapTbl.REDu[[5]]!="NO")

library(UpSetR)
list.U937_U937PMA <- list(sig.LE_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="lengthened")$gene_symbol,
                          sig.LE_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="lengthened")$gene_symbol,
                          nonsig.LE_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]>0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.LE_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]>0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol,
                          sig.SH_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[3]]=="shortened")$gene_symbol,
                          sig.SH_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[5]]=="shortened")$gene_symbol,
                          nonsig.SH_U937 = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[2]]<0 & heatmapTbl.REDu2[[3]]=="NO")$gene_symbol,
                          nonsig.SH_U937PMA = subset(heatmapTbl.REDu2,heatmapTbl.REDu2[[4]]<0 & heatmapTbl.REDu2[[5]]=="NO")$gene_symbol)
upset(fromList(list.U937_U937PMA), nsets = 8, nintersects = NA, sets = names(list.U937_U937PMA), keep.order = T, order.by = "freq", text.scale = 1.5)





List.REDi.8h3 <- lapply(list.REDi2, function(x) { x <- x[c(1,9,14)]; x })
heatmapTbl.REDi <- Reduce(function(x, y) merge(x, y), List.REDi.8h3)
heatmapTbl.REDi2 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO" | heatmapTbl.REDi[[5]]!="NO")

heatmapTbl.REDi2.U937 <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[3]]!="NO")
heatmapTbl.REDi2.U937PMA <- subset(heatmapTbl.REDi, heatmapTbl.REDi[[5]]!="NO")

library(UpSetR)
list.U937_U937PMA <- list(sig.SU_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="lengthened")$gene_symbol,
                          sig.SU_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="lengthened")$gene_symbol,
                          nonsig.SU_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]>0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.SU_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]>0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol,
                          sig.AC_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[3]]=="shortened")$gene_symbol,
                          sig.AC_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[5]]=="shortened")$gene_symbol,
                          nonsig.AC_U937 = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[2]]<0 & heatmapTbl.REDi2[[3]]=="NO")$gene_symbol,
                          nonsig.AC_U937PMA = subset(heatmapTbl.REDi2,heatmapTbl.REDi2[[4]]<0 & heatmapTbl.REDi2[[5]]=="NO")$gene_symbol)
upset(fromList(list.U937_U937PMA), nsets = 8, nintersects = NA, sets = names(list.U937_U937PMA), keep.order = T, order.by = "freq", text.scale = 1.5)


