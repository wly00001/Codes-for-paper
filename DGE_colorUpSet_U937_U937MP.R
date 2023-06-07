
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


### DGE heatmap

filenames <- list.files("./_DGE_RPM/FC2.0/", pattern="*.csv", full.names=TRUE)
filenames
ldf <- lapply(filenames[c(8,10)], read.csv)
names(ldf) <- c("U937_1_7h","U937PMA_1_7h")
# ldf <- lapply(filenames[c(7,9)], read.csv)
# names(ldf) <- c("U937_1_2h","U937PMA_1_2h")

ldf.QSF <- ldf



ldf.JTE <- c(ldf.QSF)



# ldf.JTE <- lapply(ldf.JTE, function(x) { x[[10]] <- ifelse(x[[9]] < 0.05 & x[[6]]>log2(2), "UP", ifelse(x[[9]] < 0.05 & x[[6]]<(-log2(2)), "DN", "NO")); x })  #set log2FC cutoff as log2(2)


ldf2.JTE <- lapply(ldf.JTE, function(x) { x <- x[c(1,6,10)]; x })


heatmapTbl.DGE <- Reduce(function(x, y) merge(x, y), ldf2.JTE)

# heatmapTbl.DGE2 <- subset(heatmapTbl.DGE2, (heatmapTbl.DGE2[[2]]>=log2(1.2)|heatmapTbl.DGE2[[2]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[4]]>=log2(1.2)|heatmapTbl.DGE2[[4]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[6]]>=log2(1.2)|heatmapTbl.DGE2[[6]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[8]]>=log2(1.2)|heatmapTbl.DGE2[[8]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[10]]>=log2(1.2)|heatmapTbl.DGE2[[10]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[12]]>=log2(1.2)|heatmapTbl.DGE2[[12]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[14]]>=log2(1.2)|heatmapTbl.DGE2[[14]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[16]]>=log2(1.2)|heatmapTbl.DGE2[[16]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[18]]>=log2(1.2)|heatmapTbl.DGE2[[18]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[20]]>=log2(1.2)|heatmapTbl.DGE2[[20]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[22]]>=log2(1.2)|heatmapTbl.DGE2[[22]]<=log2(1/1.2)) &
#                             (heatmapTbl.DGE2[[24]]>=log2(1.2)|heatmapTbl.DGE2[[24]]<=log2(1/1.2)))

heatmapTbl.DGE2 <- subset(heatmapTbl.DGE, heatmapTbl.DGE[[3]]!="NO" | heatmapTbl.DGE[[5]]!="NO")
heatmapTbl.DGE2.U937 <- subset(heatmapTbl.DGE, heatmapTbl.DGE[[3]]!="NO")
heatmapTbl.DGE2.U937PMA <- subset(heatmapTbl.DGE, heatmapTbl.DGE[[5]]!="NO")

library(UpSetR)
list.U937_U937PMA <- list(sig.UP_U937 = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[3]]=="UP")$gene_symbol,
                          sig.UP_U937PMA = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[5]]=="UP")$gene_symbol,
                          nonsig.UP_U937 = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[2]]>0 & heatmapTbl.DGE2[[3]]=="NO")$gene_symbol,
                          nonsig.UP_U937PMA = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[4]]>0 & heatmapTbl.DGE2[[5]]=="NO")$gene_symbol,
                          sig.DN_U937 = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[3]]=="DN")$gene_symbol,
                          sig.DN_U937PMA = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[5]]=="DN")$gene_symbol,
                          nonsig.DN_U937 = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[2]]<0 & heatmapTbl.DGE2[[3]]=="NO")$gene_symbol,
                          nonsig.DN_U937PMA = subset(heatmapTbl.DGE2,heatmapTbl.DGE2[[4]]<0 & heatmapTbl.DGE2[[5]]=="NO")$gene_symbol)
upset(fromList(list.U937_U937PMA), nsets = 8, nintersects = NA, sets = names(list.U937_U937PMA), keep.order = T, order.by = "freq", text.scale = 1.5)

# write.csv(heatmapTbl.DGE2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig7f.csv"),row.names = T)



DF <- heatmapTbl.DGE

DF$regu.adj <- "NO"
DF$regu.adj[DF$regu.adj.U937_JTE607_6h.v.U937_Mock_6h!="NO" |  DF$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h!="NO"] <- "sig"
DF$regu.adj[DF$regu.adj.U937_JTE607_6h.v.U937_Mock_6h=="UP" &  DF$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h=="UP"] <- "UP"
DF$regu.adj[DF$regu.adj.U937_JTE607_6h.v.U937_Mock_6h=="DN" &  DF$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h=="DN"] <- "DN"


names(DF)

# LPS
cor_all <- cor.test(DF$log2.U937_JTE607_6h.v.U937_Mock_6h,DF$log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h)
cor_sig <- cor.test(subset(DF, regu.adj!="NO")$log2.U937_JTE607_6h.v.U937_Mock_6h,
                    subset(DF, regu.adj!="NO")$log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h)
ggplot() +
  geom_point(data=DF, aes(x=log2.U937_JTE607_6h.v.U937_Mock_6h, y=log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h), color = "gray", size = 1) +
  geom_point(data=subset(DF, regu.adj == "sig"), aes(x=log2.U937_JTE607_6h.v.U937_Mock_6h, y=log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h), color = "black", size = 1) +
  geom_point(data=subset(DF, regu.adj == "UP"), aes(x=log2.U937_JTE607_6h.v.U937_Mock_6h, y=log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h), color = "red") +
  geom_point(data=subset(DF, regu.adj == "DN"), aes(x=log2.U937_JTE607_6h.v.U937_Mock_6h, y=log2.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h), color = "blue") +
  coord_fixed(ratio = 1)+
  theme_classic()+
  theme(axis.text.y=element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text.x=element_text(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  scale_x_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
  scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
  geom_vline(xintercept=0, col="black", linetype="solid") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  labs(title = paste("blue=",table(DF$regu.adj)[1],
                     ", red=",table(DF$regu.adj)[4],
                     ", black=",table(DF$regu.adj)[3],
                     ", gray=",table(DF$regu.adj)[2],
                     "\ncor.all=",round(cor_all$estimate,2),", pval=",formatC(cor_all$p.value, format = "e", digits = 2),
                     "\ncor.blueredblack=",round(cor_sig$estimate,2),", pval=",formatC(cor_sig$p.value, format = "e", digits = 2),
                     sep = ""), 
       x = "Log2ratio, U937_JTE607_6h.v.U937_Mock_6h", y = "Log2ratio, U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h", color = "")




















library(tidyverse)

### DGE heatmap

filenames <- list.files("./_DGE_RPM/FC2.0/", pattern="*.csv", full.names=TRUE)
filenames
ldf <- lapply(filenames[c(7:10)], read.csv)
names(ldf) <- c("U937_1_2h","U937_1_7h","U937PMA_1_2h","U937PMA_1_7h")

ldf.QSF <- ldf

b1 <- ldf.QSF[["U937_1_2h"]]
b2 <- ldf.QSF[["U937_1_7h"]]
library(UpSetR)
UPlist <- list(U937_1_2h_UP = b1[b1$regu.adj.U937_JTE607_1h.v.U937_Mock_1h=="UP",]$gene_symbol, 
               U937_1_2h_NC = b1[b1$regu.adj.U937_JTE607_1h.v.U937_Mock_1h=="NO",]$gene_symbol,
               U937_1_2h_DN = b1[b1$regu.adj.U937_JTE607_1h.v.U937_Mock_1h=="DN",]$gene_symbol,
               U937_1_7h_UP = b2[b2$regu.adj.U937_JTE607_6h.v.U937_Mock_6h=="UP",]$gene_symbol,
               U937_1_7h_NC = b2[b2$regu.adj.U937_JTE607_6h.v.U937_Mock_6h=="NO",]$gene_symbol,
               U937_1_7h_DN = b2[b2$regu.adj.U937_JTE607_6h.v.U937_Mock_6h=="DN",]$gene_symbol)
upset(fromList(UPlist), nsets = 6, nintersects = NA, sets = names(UPlist)[1:6], keep.order = T, order.by = "freq", text.scale = 1.5)

b1 <- ldf.QSF[["U937PMA_1_2h"]]
b2 <- ldf.QSF[["U937PMA_1_7h"]]
library(UpSetR)
UPlist <- list(U937MP_1_2h_UP = b1[b1$regu.adj.U937_PMA_JTE607_1h.v.U937_PMA_Mock_1h=="UP",]$gene_symbol, 
               U937MP_1_2h_NC = b1[b1$regu.adj.U937_PMA_JTE607_1h.v.U937_PMA_Mock_1h=="NO",]$gene_symbol,
               U937MP_1_2h_DN = b1[b1$regu.adj.U937_PMA_JTE607_1h.v.U937_PMA_Mock_1h=="DN",]$gene_symbol,
               U937MP_1_7h_UP = b2[b2$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h=="UP",]$gene_symbol,
               U937MP_1_7h_NC = b2[b2$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h=="NO",]$gene_symbol,
               U937MP_1_7h_DN = b2[b2$regu.adj.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h=="DN",]$gene_symbol)
upset(fromList(UPlist), nsets = 6, nintersects = NA, sets = names(UPlist)[1:6], keep.order = T, order.by = "freq", text.scale = 1.5)


