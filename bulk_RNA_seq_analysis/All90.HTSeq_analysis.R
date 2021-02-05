source("epidish.R")


ref.m <- as.matrix(read.csv("LM22_cibersort_immune_signature.txt", header=T, row.names=1, sep="\t", check.names=F))
data.m <- as.matrix(read.csv("All90.HTSeq.rpkm.xls", header=T, row.names=1, sep="\t", check.names=F))

maxit <- 50
nu.v <- c(0.25,0.5,0.75)

celltype.o <- DoCBS(data.m, ref.m, nu.v);

cellFrac.m <- t(t(celltype.o$estF))


cellFracCom.m <- matrix(NA, ncol=10, nrow=nrow(cellFrac.m))
colnames(cellFracCom.m) <- c("B_cells", "CD8T_cells", "CD4T_cells", "NK_cell", "Monocytes", "Neutrophils", "Eosinophils", "Dendritic_cells", "T_cells", "Macrophages")
rownames(cellFracCom.m) <- rownames(cellFrac.m )

bcells.idx <- c(1,2)
cd8tcells.idx <- c(4)
cd4tcells.idx <- c(5,6,7)
nkcell.idx <- c(11,12)
Monocytes.idx <- c(13)
Neutrophils.idx <- c(22)
Eosinophils.idx <- c(21)
Dendritic.idx <- c(17,18)
Tcells.idx <- 4:10
Macrophages.idx <- 14:16

cellFracCom.m[,1] <- rowSums(cellFrac.m[,bcells.idx])
cellFracCom.m[,2] <- cellFrac.m[,cd8tcells.idx]
cellFracCom.m[,3] <- rowSums(cellFrac.m[,cd4tcells.idx])
cellFracCom.m[,4] <- rowSums(cellFrac.m[,nkcell.idx])
cellFracCom.m[,5] <- cellFrac.m[,Monocytes.idx]
cellFracCom.m[,6] <- cellFrac.m[,Neutrophils.idx]
cellFracCom.m[,7] <- cellFrac.m[,Eosinophils.idx]
cellFracCom.m[,8] <- rowSums(cellFrac.m[,Dendritic.idx])
cellFracCom.m[,9] <- rowSums(cellFrac.m[,Tcells.idx])
cellFracCom.m[,10] <- rowSums(cellFrac.m[,Macrophages.idx])

cellFracCom.m <- cellFracCom.m * 100


write.table(cellFrac.m, file="DATA_cellFrac_Original.txt", sep="\t", quote=F, col.names=NA)

write.table(cellFracCom.m, file="DATA_cellFracCom.txt", sep="\t", quote=F, col.names=NA)


PhenoTypes.df <- read.csv("SampleClass.txt", header=F, sep="\t")
tmp.idx <- match(rownames(cellFracCom.m), PhenoTypes.df[,1])
PhenoTypes.df <- PhenoTypes.df[tmp.idx,]

cancer.idx <- which(PhenoTypes.df$V2 == "cancer")
cellFracCom.m <- cellFracCom.m[cancer.idx,]
data.m <- data.m[,cancer.idx]
############

cellFracCom.df <- as.data.frame(cellFracCom.m)
cellFracCom.df$Index <- NA

BCL9.idx <- which(rownames(data.m) == "BCL9")
BCL9.v <- data.m[BCL9.idx,]
#fivenum(BCL9.v)
quantile(BCL9.v)
cellFracCom.df$Index[which(BCL9.v <= quantile(BCL9.v)[2] )] <- "Lower"
cellFracCom.df$Index[which(BCL9.v >= quantile(BCL9.v)[4] )] <- "Higher"

library(reshape2)
library(ggplot2)
library(ggpubr)
data.df <- melt(cellFracCom.df)
data.df <- data.df[which(!is.na(data.df[,1])),]
ggplot(data.df, aes(x = Index, y=value)) + geom_violin() + facet_wrap(variable ~ . ) + theme_bw()
ggsave("tmp.pdf")
ggviolin(data = data.df ,x='Index',y='value', fill='Index', palette = c("#FC4E07", "#00AFBB"), facet.by = 'variable',add = "boxplot",add.params = list(fill="white"), ggtheme = theme_light(), legend = "", ylab="Cibersort infiltration estimates", xlab="")+ stat_compare_means(comparisons=list(c('Lower','Higher')),label = "p.signif")
ggsave("DATA_results/DATA_BCL9_Cibersort_Plot.pdf")

#############
cellFracCom.df <- as.data.frame(cellFracCom.m)
cellFracCom.df$Index <- NA

BCL9L.idx <- which(rownames(data.m) == "BCL9L")
BCL9L.v <- data.m[BCL9L.idx,]
#fivenum(BCL9L.v)
quantile(BCL9L.v)
cellFracCom.df$Index[which(BCL9L.v <= quantile(BCL9L.v)[2] )] <- "Lower"
cellFracCom.df$Index[which(BCL9L.v >= quantile(BCL9L.v)[4] )] <- "Higher"


data.df <- melt(cellFracCom.df)
data.df <- data.df[which(!is.na(data.df[,1])),]
#ggplot(data.df, aes(x = Index, y=value)) + geom_violin() + facet_wrap(variable ~ . ) + theme_bw()
#ggsave("tmp.pdf")
ggviolin(data = data.df ,x='Index',y='value', fill='Index', palette = c("#FC4E07", "#00AFBB"), facet.by = 'variable',add = "boxplot",add.params = list(fill="white"), ggtheme = theme_light(), legend = "", ylab="Cibersort infiltration estimates", xlab="")+ stat_compare_means(comparisons=list(c('Lower','Higher')),label = "p.signif")
ggsave("DATA_results/DATA_BCL9L_Cibersort_Plot.pdf")

##################################################################################################

#allmutSample.v <- substr(as.vector(read.csv("Allsamples.txt", header=F,)[,1]), 1, 16)
#CTNNB1_Sample.v <- substr(as.vector(read.csv("CTNNB1_mutated_samples.txt", header=F,)[,1]), 1, 16)

CTNNB1_Sample.v <- c("A1076", "A1099", "LV0042T", "LV0046T", "LV0068T", "LV0069T", "LV0092T", "LV61101T", "T15_422", "T16", "T275")

cellFracCom.df <- as.data.frame(cellFracCom.m)
cellFracCom.df$Index <- "WT"


tmp.idx <- match(rownames(cellFracCom.df), CTNNB1_Sample.v)
cellFracCom.df[which(!is.na(tmp.idx)), 10] <- "MUT"


data.df <- melt(cellFracCom.df)
data.df <- data.df[which(!is.na(data.df[,1])),]
#ggplot(data.df, aes(x = Index, y=value)) + geom_violin() + facet_wrap(variable ~ . ) + theme_bw()
#ggsave("tmp.pdf")
ggviolin(data = data.df ,x='Index',y='value', fill='Index', palette = c("#E7B800", "#3366CC"), facet.by = 'variable',add = "boxplot",add.params = list(fill="white"), ggtheme = theme_light(), legend = "", ylab="Cibersort infiltration estimates", xlab="")+ stat_compare_means(comparisons=list(c('MUT', 'WT')),label = "p.signif")
ggsave("DATA_results/DATA_CTNNB1_mut_Cibersort_Plot.pdf")





