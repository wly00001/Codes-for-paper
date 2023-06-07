
workplace <- "//cifs.wistar.upenn.edu/Tian"
# workplace <- "/run/user/1000/gvfs/smb-share:server=cifs.wistar.upenn.edu,share=tian"


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




filenames <- list.files(paste0(workplace,"/linux/lwang/drive/project/21029-09/_DGE/FC1.2/"), pattern="*.csv", full.names=TRUE)
filenames

names(filenames) <- gsub("YC_103_pA_1_1","24DnJn",names(filenames))
names(filenames) <- gsub("YC_103_pA_1_2","24DnJp",names(filenames))
names(filenames) <- gsub("YC_103_pA_1_3","24DpJn",names(filenames))
names(filenames) <- gsub("YC_103_pA_1_4","24DpJp",names(filenames))
names(filenames) <- gsub("YC_103_pA_2_5","48DnJn",names(filenames))
names(filenames) <- gsub("YC_103_pA_2_6","48DnJp",names(filenames))
names(filenames) <- gsub("YC_103_pA_2_7","48DpJn",names(filenames))
names(filenames) <- gsub("YC_103_pA_2_8","48DpJp",names(filenames))

GEplot_MA <- list()
GEplot_histogram <- list()

for (k in 1:length(filenames)) {
  
  # k=7
  
  treat_gene_expression5 <- read.table(filenames[k], header = T, sep = ",")[c(1,4:6,9,10)]
  treat_gene_expression5$meanRPM <- (treat_gene_expression5[[2]]+treat_gene_expression5[[3]])/2
  treat_gene_expression5[c(2,3)] <-NULL
  names(treat_gene_expression5) <- c("gene_symbol","L2FC","pval.fisher.adj","regu.adj","meanRPM")
  treat_gene_expression5 <- treat_gene_expression5 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
  treat_gene_expression5 <- na.omit(treat_gene_expression5)
  
  # treat_gene_expression5[["regu.adj"]] <- ifelse(treat_gene_expression5[[paste0("pval.fisher.adj")]] < 0.05 & treat_gene_expression5[[paste0("L2FC")]]>log2(2), "UP", ifelse(treat_gene_expression5[[paste0("pval.fisher.adj")]] < 0.05 & treat_gene_expression5[[paste0("L2FC")]]<(-log2(2)), "DN", "NO")) #set log2FC cutoff as log2(2)
  
  sample <- gsub(".BH.csv","",tail(strsplit(filenames[k],"/")[[1]],1))
  sample <- gsub("YC_103_pA_1_1","24DnJn",sample)
  sample <- gsub("YC_103_pA_1_2","24DnJp",sample)
  sample <- gsub("YC_103_pA_1_3","24DpJn",sample)
  sample <- gsub("YC_103_pA_1_4","24DpJp",sample)
  sample <- gsub("YC_103_pA_2_5","48DnJn",sample)
  sample <- gsub("YC_103_pA_2_6","48DnJp",sample)
  sample <- gsub("YC_103_pA_2_7","48DpJn",sample)
  sample <- gsub("YC_103_pA_2_8","48DpJp",sample)
  
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
    # scale_x_continuous(limits = c(0, 15), breaks = c(0,5,10,15))+
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
    geom_hline(yintercept=c(-1,-log2(1.2),log2(1.2),1), col="black",linetype="dashed") +
    ggtitle(paste("DGE, ",sample,
                  "\nDN=",table(treat_gene_expression5$regu.adj)[1],
                  ", NC=",table(treat_gene_expression5$regu.adj)[2],
                  ", UP=",table(treat_gene_expression5$regu.adj)[3],
                  ", UP/DN=",round(table(treat_gene_expression5$regu.adj)[3]/table(treat_gene_expression5$regu.adj)[1],2),sep = ""))
  
  
  DF_histogram <- subset(treat_gene_expression5, regu.adj != "NO")
  GEplot_histogram[[sample]] <- ggplot(DF_histogram, aes(x=log2(meanRPM))) + 
    geom_histogram(color="black", fill="white", size=1.2) +
    theme_classic()+
    scale_x_continuous(limits = c(0, 10), breaks = c(0,5,10))+
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
    ggtitle(paste("# of regulated genes, ",sample, sep = ""))
  
  
  # write.csv(treat_gene_expression5, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig9f.csv"),row.names = T)
  
  
}



library(cowplot)
plot_grid(plotlist=GEplot_MA[c(1,5,6,10)],ncol = 2)
plot_grid(plotlist=GEplot_MA[c(2,7)],ncol = 1)

# plot_grid(plotlist=c(GEplot_histogram),ncol = 5)



library(tidyverse)

### DGE heatmap

filenames <- list.files(paste0(workplace,"/linux/lwang/drive/project/21029-09/_DGE/FC1.2/"), pattern="*.csv", full.names=TRUE)
filenames
ldf <- lapply(filenames, read.csv)
names(ldf) <- substr(filenames,80,nchar(filenames)-7)
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


b1 <- new.ldf[["48DnJp.v.48DnJn"]]
b2 <- new.ldf[["48DpJp.v.48DpJn"]]
library(UpSetR)
UPlist <- list(DGE_48DnJp.v.48DnJn_UP = b1[b1$regu.adj.48DnJp.v.48DnJn=="UP",]$gene_symbol, 
               DGE_48DnJp.v.48DnJn_NC = b1[b1$regu.adj.48DnJp.v.48DnJn=="NO",]$gene_symbol,
               DGE_48DnJp.v.48DnJn_DN = b1[b1$regu.adj.48DnJp.v.48DnJn=="DN",]$gene_symbol,
               DGE_48DpJp.v.48DpJn_UP = b2[b2$regu.adj.48DpJp.v.48DpJn=="UP",]$gene_symbol,
               DGE_48DpJp.v.48DpJn_NC = b2[b2$regu.adj.48DpJp.v.48DpJn=="NO",]$gene_symbol,
               DGE_48DpJp.v.48DpJn_DN = b2[b2$regu.adj.48DpJp.v.48DpJn=="DN",]$gene_symbol)
upset(fromList(UPlist), nsets = 6, nintersects = NA, sets = names(UPlist)[1:6], keep.order = T, order.by = "freq", text.scale = 1.5)

# write.csv(merge(new.ldf[["48DnJp.v.48DnJn"]][c(1,10)],new.ldf[["48DpJp.v.48DpJn"]][c(1,10)]), paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig9i.csv"),row.names = T)


# comparing Log2Ratio between U937 and U937+PMA

##1h
DF <- merge(b1,b2,by="gene_symbol")
names(DF)


DF$regu.adj <- "NO"
DF$regu.adj[DF$regu.adj.48DnJp.v.48DnJn!="NO" |  DF$regu.adj.48DpJp.v.48DpJn!="NO"] <- "sig"
DF$regu.adj[DF$regu.adj.48DnJp.v.48DnJn=="UP" &  DF$regu.adj.48DpJp.v.48DpJn=="UP"] <- "UP"
DF$regu.adj[DF$regu.adj.48DnJp.v.48DnJn=="DN" &  DF$regu.adj.48DpJp.v.48DpJn=="DN"] <- "DN"


names(DF)

# LPS
cor_all <- cor.test(DF$log2.48DnJp.v.48DnJn,DF$log2.48DpJp.v.48DpJn)
cor_sig <- cor.test(subset(DF, regu.adj!="NO")$log2.48DnJp.v.48DnJn,
                    subset(DF, regu.adj!="NO")$log2.48DpJp.v.48DpJn)
ggplot() +
  geom_point(data=DF, aes(x=log2.48DnJp.v.48DnJn, y=log2.48DpJp.v.48DpJn), color = "gray", size = 1) +
  geom_point(data=subset(DF, regu.adj == "sig"), aes(x=log2.48DnJp.v.48DnJn, y=log2.48DpJp.v.48DpJn), color = "black", size = 1) +
  geom_point(data=subset(DF, regu.adj == "UP"), aes(x=log2.48DnJp.v.48DnJn, y=log2.48DpJp.v.48DpJn), color = "red") +
  geom_point(data=subset(DF, regu.adj == "DN"), aes(x=log2.48DnJp.v.48DnJn, y=log2.48DpJp.v.48DpJn), color = "blue") +
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
       x = "Log2ratio, 48DnJp.v.48DnJn,", y = "Log2ratio, 48DpJp.v.48DpJn", color = "")

# write.csv(DF[c(1,6,10,15,19,20)], paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/fig9h.csv"),row.names = T)



b1 <- new.ldf[["24DnJp.v.24DnJn"]]
b2 <- new.ldf[["24DpJp.v.24DpJn"]]
library(UpSetR)
UPlist <- list(DGE_24DnJp.v.24DnJn_UP = b1[b1$regu.adj.24DnJp.v.24DnJn=="UP",]$gene_symbol, 
               DGE_24DnJp.v.24DnJn_NC = b1[b1$regu.adj.24DnJp.v.24DnJn=="NO",]$gene_symbol,
               DGE_24DnJp.v.24DnJn_DN = b1[b1$regu.adj.24DnJp.v.24DnJn=="DN",]$gene_symbol,
               DGE_24DpJp.v.24DpJn_UP = b2[b2$regu.adj.24DpJp.v.24DpJn=="UP",]$gene_symbol,
               DGE_24DpJp.v.24DpJn_NC = b2[b2$regu.adj.24DpJp.v.24DpJn=="NO",]$gene_symbol,
               DGE_24DpJp.v.24DpJn_DN = b2[b2$regu.adj.24DpJp.v.24DpJn=="DN",]$gene_symbol)
upset(fromList(UPlist), nsets = 6, nintersects = NA, sets = names(UPlist)[1:6], keep.order = T, order.by = "freq", text.scale = 1.5)



# comparing Log2Ratio between U937 and U937+PMA

##1h
DF <- merge(b1,b2,by="gene_symbol")
names(DF)


DF$regu.adj <- "NO"
DF$regu.adj[DF$regu.adj.24DnJp.v.24DnJn!="NO" |  DF$regu.adj.24DpJp.v.24DpJn!="NO"] <- "sig"
DF$regu.adj[DF$regu.adj.24DnJp.v.24DnJn=="UP" &  DF$regu.adj.24DpJp.v.24DpJn=="UP"] <- "UP"
DF$regu.adj[DF$regu.adj.24DnJp.v.24DnJn=="DN" &  DF$regu.adj.24DpJp.v.24DpJn=="DN"] <- "DN"


names(DF)

# LPS
cor_all <- cor.test(DF$log2.24DnJp.v.24DnJn,DF$log2.24DpJp.v.24DpJn)
cor_sig <- cor.test(subset(DF, regu.adj!="NO")$log2.24DnJp.v.24DnJn,
                    subset(DF, regu.adj!="NO")$log2.24DpJp.v.24DpJn)
ggplot() +
  geom_point(data=DF, aes(x=log2.24DnJp.v.24DnJn, y=log2.24DpJp.v.24DpJn), color = "gray", size = 1) +
  geom_point(data=subset(DF, regu.adj == "sig"), aes(x=log2.24DnJp.v.24DnJn, y=log2.24DpJp.v.24DpJn), color = "black", size = 1) +
  geom_point(data=subset(DF, regu.adj == "UP"), aes(x=log2.24DnJp.v.24DnJn, y=log2.24DpJp.v.24DpJn), color = "red") +
  geom_point(data=subset(DF, regu.adj == "DN"), aes(x=log2.24DnJp.v.24DnJn, y=log2.24DpJp.v.24DpJn), color = "blue") +
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
       x = "Log2ratio, 24DnJp.v.24DnJn,", y = "Log2ratio, 24DpJp.v.24DpJn", color = "")




