
# 700*700 .SVG

library(DESeq2)
library(ggplot2)
library(reshape)
library(tidyverse)
library(ggrepel)



REDu_NP <- read.table(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/APA_JTE_PCF11kd_FIP1oe/3READSpy APA analysis step2 - replicated - TPA.tbl"), header = T, sep = "\t")[c(1,44)]
names(REDu_NP)[c(2)] <- "REDu_NP"
REDi_NP <- read.table(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/APA_JTE_PCF11kd_FIP1oe/3READSpy APA analysis step2 - replicated - IPA.tbl"), header = T, sep = "\t")[c(1,44)]
names(REDi_NP)[c(2)] <- "REDi_NP"


List.REDu <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/APA_JTE_PCF11kd_FIP1oe/List.REDu.RData"))
List.REDi <- readRDS(paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/CPA factors/PCF11/Kamieniarz-Gdula_2019_Molecular_Cell/APA_JTE_PCF11kd_FIP1oe/List.REDi.RData"))

List.REDu2 <- list()
List.REDi2 <- list()

for (sample in names(List.REDu)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDu[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDu2[[sample]] <- df
  
}

for (sample in names(List.REDi)) {
  
  # sample="HeLa_JTE_01"
  
  df <- List.REDi[[sample]]
  names(df)[c(5:14)] <- paste0(names(df)[c(5:14)],".",sample)
  List.REDi2[[sample]] <- df
  
}

List.REDu2 <- lapply(List.REDu2, function(x) { x <- x[c(1,9)]; x })
List.REDi2 <- lapply(List.REDi2, function(x) { x <- x[c(1,9)]; x })

List.REDu3 <- c(List.REDu2,list(REDu_NP))
List.REDi3 <- c(List.REDi2,list(REDi_NP))

List.REDu3 <- lapply(List.REDu3, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu3 <- lapply(List.REDu3, function(x) { x<-na.omit(x); x })
List.REDi3 <- lapply(List.REDi3, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDi3 <- lapply(List.REDi3, function(x) { x<-na.omit(x); x })


DF_REDu <- Reduce(function(x, y) merge(x, y), List.REDu3[-c(7,9,10)])
DF_REDi <- Reduce(function(x, y) merge(x, y), List.REDi3[-c(7,9,10)])






tbl <- data.frame(DF_REDu[-c(1)],row.names=DF_REDu$gene_symbol)
tbl <- tbl %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
tbl <- na.omit(tbl)

tbl2 <- tbl[c(1:4,14)]
# write.csv(tbl2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS4a.csv"),row.names = T)

x <- c("REDu1.HeLa_JTE_1.v.HeLa_DMSO","REDu1.HeLa_JTE_1.v.HeLa_DMSO","REDu1.HeLa_JTE_1.v.HeLa_DMSO","REDu1.HeLa_JTE_1.v.HeLa_DMSO",
       "REDu1.HeLa_JTE_10.v.HeLa_DMSO","REDu1.HeLa_JTE_10.v.HeLa_DMSO","REDu1.HeLa_JTE_10.v.HeLa_DMSO",
       "REDu1.HepG2_JTE1.v.HepG2_DMSO","REDu1.HepG2_JTE1.v.HepG2_DMSO",
       "REDu1.HepG2_JTE10.v.HepG2_DMSO")
y <- c("REDu1.HeLa_JTE_10.v.HeLa_DMSO","REDu1.HepG2_JTE1.v.HepG2_DMSO","REDu1.HepG2_JTE10.v.HepG2_DMSO","REDu_NP",
       "REDu1.HepG2_JTE1.v.HepG2_DMSO","REDu1.HepG2_JTE10.v.HepG2_DMSO","REDu_NP",
       "REDu1.HepG2_JTE10.v.HepG2_DMSO","REDu_NP",
       "REDu_NP")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 50) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)

tbl2 <- tbl[c(3:6)]
# write.csv(tbl2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS6c.csv"),row.names = T)

x <- c("REDu1.HepG2_JTE1.v.HepG2_DMSO","REDu1.HepG2_JTE1.v.HepG2_DMSO","REDu1.HepG2_JTE1.v.HepG2_DMSO",
       "REDu1.HepG2_JTE10.v.HepG2_DMSO","REDu1.HepG2_JTE10.v.HepG2_DMSO",
       "REDu1.U937_JTE607_6h.v.U937_Mock_6h")
y <- c("REDu1.HepG2_JTE10.v.HepG2_DMSO","REDu1.U937_JTE607_6h.v.U937_Mock_6h","REDu1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "REDu1.U937_JTE607_6h.v.U937_Mock_6h","REDu1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "REDu1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 50) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)







tbl <- data.frame(DF_REDi[-c(1)],row.names=DF_REDi$gene_symbol)
tbl <- tbl %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
tbl <- na.omit(tbl)

tbl2 <- tbl[c(1:4,14)]
# write.csv(tbl2, paste0("//cifs.wistar.upenn.edu/Tian/linux/lwang/drive/project/21029-07/allJTE_06_07/sum0313paper/email/final/Source Data/figS4b.csv"),row.names = T)

x <- c("REDi1.HeLa_JTE_1.v.HeLa_DMSO","REDi1.HeLa_JTE_1.v.HeLa_DMSO","REDi1.HeLa_JTE_1.v.HeLa_DMSO","REDi1.HeLa_JTE_1.v.HeLa_DMSO",
       "REDi1.HeLa_JTE_10.v.HeLa_DMSO","REDi1.HeLa_JTE_10.v.HeLa_DMSO","REDi1.HeLa_JTE_10.v.HeLa_DMSO",
       "REDi1.HepG2_JTE1.v.HepG2_DMSO","REDi1.HepG2_JTE1.v.HepG2_DMSO",
       "REDi1.HepG2_JTE10.v.HepG2_DMSO")
y <- c("REDi1.HeLa_JTE_10.v.HeLa_DMSO","REDi1.HepG2_JTE1.v.HepG2_DMSO","REDi1.HepG2_JTE10.v.HepG2_DMSO","REDi_NP",
       "REDi1.HepG2_JTE1.v.HepG2_DMSO","REDi1.HepG2_JTE10.v.HepG2_DMSO","REDi_NP",
       "REDi1.HepG2_JTE10.v.HepG2_DMSO","REDi_NP",
       "REDi_NP")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 50) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)

tbl2 <- tbl[c(3:6)]
x <- c("REDi1.HepG2_JTE1.v.HepG2_DMSO","REDi1.HepG2_JTE1.v.HepG2_DMSO","REDi1.HepG2_JTE1.v.HepG2_DMSO",
       "REDi1.HepG2_JTE10.v.HepG2_DMSO","REDi1.HepG2_JTE10.v.HepG2_DMSO",
       "REDi1.U937_JTE607_6h.v.U937_Mock_6h")
y <- c("REDi1.HepG2_JTE10.v.HepG2_DMSO","REDi1.U937_JTE607_6h.v.U937_Mock_6h","REDi1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "REDi1.U937_JTE607_6h.v.U937_Mock_6h","REDi1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h",
       "REDi1.U937_PMA_JTE607_6h.v.U937_PMA_Mock_6h")

plot <- list()

for (n in 1:length(x)){
  
  # n=1
  
  Nx <- x[n]
  Ny <- y[n]
  
  
  plot[[n]] <- ggplot(data=tbl2, aes(x=.data[[Nx]], y=.data[[Ny]])) +
    geom_point() +
    geom_bin2d(bins = 50) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    geom_vline(xintercept=0, col="black", linetype="solid") +
    geom_hline(yintercept=0, col="black", linetype="solid") +
    labs(x = Nx, y = Ny)
  
}

library(cowplot)
plot_grid(plotlist=plot,ncol = 2)


















