
library(ggpubr)

clust0.idx <- which(Idents(sce) == 0)

data1.m <- data.m[,clust0.idx]

metadata1.df <- metadata.df[clust0.idx, ]

table(metadata1.df[,2],metadata1.df[,3])

tmp2.m <- rbind(as.vector(metadata1.df[,1], data1.m)
###########
tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-PD1"), grep("Pre", metadata1.df[,1])) , which(data1.m[7840, ] >0) )


CD27exp.v <- as.vector(data1.m[7840,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])

tmp.df <- as.data.frame(cbind(CD27exp.v, class.v))
tmp.df[,1] <- as.numeric(CD27exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("CD27_exp_boxplot_by_PD1.pdf", width=5, height=5)

write.table(tmp.df , file="CD27_exp_boxplot_by_PD1.csv", quote=F, sep="\t", row.names=F)



tmp.df$CD27class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$CD27class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- tmp1.df/colSums(tmp1.df)

pdf("CD27_exp_barplot_by_PD1.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()


###########
tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-CTLA4"), grep("Pre", metadata1.df[,1])) , which(data1.m[7840, ] >0) )

CD27exp.v <- as.vector(data1.m[7840,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])
#boxplot(CD27exp.v ~ class.v)
#t.test(CD27exp.v ~ class.v)
tmp.df <- as.data.frame(cbind(CD27exp.v, class.v))
tmp.df[,1] <- as.numeric(CD27exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("CD27_exp_boxplot_by_CTLA4.pdf", width=5, height=5)


tmp.df$CD27class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$CD27class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- tmp1.df/colSums(tmp1.df)

pdf("CD27_exp_barplot_by_CTLA4.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()


###########
tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-CTLA4+PD1"), grep("Pre", metadata1.df[,1])) , which(data1.m[7840, ] >0) )

CD27exp.v <- as.vector(data1.m[7840,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])
#boxplot(CD27exp.v ~ class.v)
#t.test(CD27exp.v ~ class.v)
tmp.df <- as.data.frame(cbind(CD27exp.v, class.v))
tmp.df[,1] <- as.numeric(CD27exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("CD27_exp_boxplot_by_CTLA4_PD1.pdf", width=5, height=5)

tmp.df$CD27class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$CD27class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- tmp1.df/colSums(tmp1.df)


pdf("CD27_exp_barplot_by_CTLA4_PD1.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()


#################
#################
 


###########
#tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-PD1"), grep("Pre", metadata1.df[,1])) , which(data1.m[4562, ] >0) )
tmp1.idx <- intersect(which(metadata1.df[,3] == "anti-PD1"),  which(data1.m[4562, ] >0) )
BCL9exp.v <- as.vector(data1.m[4562,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])

tmp.df <- as.data.frame(cbind(BCL9exp.v, class.v))
tmp.df[,1] <- as.numeric(BCL9exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("BCL9_exp_boxplot_by_PD1.pdf", width=5, height=5)

write.table(tmp.df , file="BCL9_exp_boxplot_by_PD1.csv", quote=F, sep="\t", row.names=F)



tmp.df$BCL9class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$BCL9class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- round(tmp1.df/colSums(tmp1.df), digits=3)

pdf("BCL9_exp_barplot_by_PD1.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()


###########
#tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-CTLA4"), grep("Pre", metadata1.df[,1])) , which(data1.m[4562, ] >0) )
tmp1.idx <- intersect(which(metadata1.df[,3] == "anti-CTLA4"),  which(data1.m[4562, ] >0) )
BCL9exp.v <- as.vector(data1.m[4562,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])
#boxplot(BCL9exp.v ~ class.v)
#t.test(BCL9exp.v ~ class.v)
tmp.df <- as.data.frame(cbind(BCL9exp.v, class.v))
tmp.df[,1] <- as.numeric(BCL9exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("BCL9_exp_boxplot_by_CTLA4.pdf", width=5, height=5)


tmp.df$BCL9class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$BCL9class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- tmp1.df/colSums(tmp1.df)

pdf("BCL9_exp_barplot_by_CTLA4.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()


###########
#tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-CTLA4+PD1"), grep("Pre", metadata1.df[,1])) , which(data1.m[4562, ] >0) )
tmp1.idx <- intersect(which(metadata1.df[,3] == "anti-CTLA4+PD1"),  which(data1.m[4562, ] >0) )

BCL9exp.v <- as.vector(data1.m[4562,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])
#boxplot(BCL9exp.v ~ class.v)
#t.test(BCL9exp.v ~ class.v)
tmp.df <- as.data.frame(cbind(BCL9exp.v, class.v))
tmp.df[,1] <- as.numeric(BCL9exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("BCL9_exp_boxplot_by_CTLA4_PD1.pdf", width=5, height=5)

tmp.df$BCL9class <- ifelse(tmp.df$Expression > median(tmp.df$Expression), "high", "low")
tmp1.df <- table( tmp.df$Group, tmp.df$BCL9class)
p.v <- format(fisher.test(tmp1.df)$p.v, digits=3)
tmp1.df <- tmp1.df/colSums(tmp1.df)


pdf("BCL9_exp_barplot_by_CTLA4_PD1.pdf", width=5, height=5)
par(mar=c(4,4,3,9))
#par(mai=c(0.5,0.5,0.5,2))
barplot(tmp1.df[,2:1],col=c("darkolivegreen", "darkgoldenrod"), width = 0.8, space =0.8, )
legend(x=3.0,y=0.8, legend=c("Non-responder", "Responder"), fill=c("darkolivegreen", "darkgoldenrod"), xpd=TRUE  )
text(x=3.7,y=0.5, label=paste0("Fisher P=", p.v), xpd=T)
dev.off()




#################
#################
 


###########

tmp1.idx <- intersect(intersect(which(metadata1.df[,3] == "anti-PD1"), grep("Pre", metadata1.df[,1])) , which(data1.m[3886, ] >0) )

FOXM1exp.v <- as.vector(data1.m[3886,tmp1.idx])
class.v <- as.vector(metadata1.df[tmp1.idx,2])

tmp.df <- as.data.frame(cbind(FOXM1exp.v, class.v))
tmp.df[,1] <- as.numeric(FOXM1exp.v)
colnames(tmp.df) <- c("Expression", "Group")

p <- ggboxplot(tmp.df, x="Group", y="Expression" , fill="Group",  palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab="",  legend = "null" ) + stat_compare_means(method = "t.test")
p
ggsave("FOXM1_exp_boxplot_by_PD1.pdf", width=5, height=5)



write.table(tmp.df , file="FOXM1_exp_boxplot_by_PD1.csv", quote=F, sep="\t", row.names=F)




#################
#################
 
sce1 <- CreateSeuratObject(counts = data1.m, assay = "RNA", meta.data=metadata1.df, min.cells=0, min.features=0,  project = "Melanoma") 


sce1 <- FindVariableFeatures(sce1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sce1)
sce1 <- ScaleData(sce1, features = all.genes)

### 线性降维
#在scale.data上运行PCA，默认使用的基因为高变异基因
sce1 <- RunPCA(sce1, features = VariableFeatures(object = sce1))









