library(dplyr)
library(magrittr)
library(clusterProfiler)
library(enrichplot)

method_stat = "pval"

file_input     = "input.txt"
file_term2gene = "term2gene.txt"
#file_term2name = "term2name.txt"

file_out       = "gsea.xls"
file_outdir    = "gsea_plot"

geneList <- read.table( file_input, header = TRUE, sep="\t", comment.char = "", check.names=FALSE, quote = "")
colnames(geneList) <- c("gene","Fold.Change")
geneList.sort <- arrange(geneList, desc(Fold.Change));
gene <- geneList.sort$gene

glist <- geneList[,2]
names(glist) <- as.character(geneList[,1])
glist <- sort(glist,decreasing = T)

t2g <- read.table( file_term2gene, header = TRUE, sep="\t", comment.char = "", check.names=FALSE, quote = "")
#t2n <- read.table( file_term2name, header = TRUE, sep="\t", comment.char = "", check.names=FALSE, quote = "")

#enrich <- enricher(gene, TERM2GENE=t2g, TERM2NAME = t2n)

#gsea <- GSEA(glist, TERM2GENE=t2g, TERM2NAME=t2n, verbose=T, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH")
gsea <- GSEA(glist, TERM2GENE=t2g, verbose=T, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH")

res  <- gsea[,-8]
colnames(res)[1] = "NAME"
colnames(res)[6] = "P_value"
colnames(res)[7] = "P_adjust"
write.table(res, file=file_out, quote=F, col.names=T, row.names=F,sep="\t")

for(i in 1:nrow(gsea)){
  plot_flag = 0
  if(method_stat == "pval" ){
    if(gsea$pvalue[i] <= 0.05){
      p <- gseaplot2(gsea, geneSetID = i, title = gsea$Description[i], pvalue_table=T)
      plot_flag = 1
    }
  }else{
    if(gsea$p.adjust[i] <= 0.05){
      p <- gseaplot2(gsea, geneSetID = i, title = gsea$Description[i], pvalue_table=T)
      plot_flag = 1
    }
  }
  if(i > 50){
	break;
  }
  ii = sprintf("%04d", i)
  if(plot_flag >0){
    pdf(paste0(file_outdir,"/",ii,".",gsea$ID[i],".pdf"),width=8,height=8)
    print(p)
    dev.off()
  }
}

sig.idx <- which(res[,6] < 0.05)
res.sig <- res[sig.idx,]

up.idx <- which(res.sig[,5] > 0)
res.up <- res.sig[up.idx,]

down.idx <- which(res.sig[,5] < 0)
res.down <- res.sig[down.idx,]

write.table(res.up, file="gsea_up.xls", quote=F, col.names=T, row.names=F,sep="\t")
write.table(res.down, file="gsea_down.xls", quote=F, col.names=T, row.names=F,sep="\t")





