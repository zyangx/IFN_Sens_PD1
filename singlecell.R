library(Seurat)
library(ggplot2)
library(configr)
library(RcppTOML)
library(logging)
library(plyr)
library(dplyr)
library(tidyr) #### 2020-0721
library(SingleR)### 2020-0727
library(RColorBrewer)### 2020-0727

options(future.globals.maxSize = 400000 * 1024^2)
basicConfig(level='FINEST')
basicConfig(level='INFO')
suppressPackageStartupMessages(library(optparse))
savepdf<- function(outpath,outf,ispdf=1,ispng=1,notgg=0,widthle=NA,heightle =NA){
  if(notgg==0){
    if (ispdf){
      
      ggsave(outf,file=paste0(outpath,".pdf",collapse = ""),device="pdf",width=widthle,height =heightle)
    }
    if(ispng){
      ggsave(outf,file=paste0(outpath,".png",sep=""),device="png",width=widthle,height =heightle)
    }
  }
  else{
    if (ispdf){
      
      savePlot(filename=paste0(outpath,".pdf",collapse = ""),type = "png",device=dev.cur())
    }
    if(ispng){
      #ggsave(outf,file=paste0(outpath,".png",sep=""),device="png",width=5,height=5)
      loginfo("pass")
    }
    
  }
  loginfo("savefig at %s",outpath)
}

isstep <-function(x){
  outb=FALSE  
  c=paste0("outb=csc$run$step",x)
  tryCatch(
    {
      eval(parse(text=c))
      if(is.null(outb)){
        outb=FALSE
      }
      
    },error = function ( e ) {
      outb=FALSE
      logerror(print(outb))
      return(outb)
      
    },finally = {
      return(outb)
    }
  )
}

getrds <-function(x){
  backobsig=1
  c=paste0(outtemp,"step",x,".RDS")
  if(file.exists(c)){
    loginfo("import %s",c)
    backob=readRDS(c)
    backobsig=0
  }
  else if(!is.null( csc$tempdata$tempdata)){
    if(file.exists(csc$tempdata$tempdata)){
      loginfo("import otherdata %s",csc$tempdata$tempdata)
      backob=readRDS(csc$tempdata$tempdata)
      backobsig=0
    }
    
    
  }
  if(backobsig==1){
    logerror("can't find file exist in %s or %s",c,csc$tempdata$tempdata)
    stop()
  }
  return(backob)
}


option_list=list(make_option(c("-i","--infile"),  action="store", help='The input file.'))
opt<-parse_args(OptionParser(usage="%prog [options] file\n",option_list=option_list))
opt$infile="/data/haplox/users/liyq/singlecell/cell_cd27/confige.ini"
#opt$infile="/data/haplox/users/liyq/singlecell/deadline/confige.ini"   shielded by zhangdx
csc=parseTOML(opt$infile)
outdir=csc$outpath$outpath
numsap=length(csc$indata)
dir.create(outdir)
dir.create(paste(outdir,"temp",sep="/"))
outtemp=paste(outdir,"temp","",sep="/")
if(isstep(1)){
  loginfo("step1 start")
  dir.create(paste(outdir,"step1",sep="/"))
  listobjall=list()
  cat("name\tbefore\tafter\tfiltered\n",file = paste(outdir,"step1","summary.xls",sep="/"),fill= F, labels=NULL, append=F)
  for (i in 1:length(csc$indata))
  {
    csc$step1$nFeature_RNA=as.numeric(csc$step1$nFeature_RNA)
    csc$step1$percent_mt=as.numeric(csc$step1$percent_mt)
    
    dir.create(paste(outdir,"step1",names(csc$indata[i]),sep="/"))
    if(csc$step1$filetype=="10x"){
    pbmc_data= Read10X(data.dir = csc$indata[[i]])
    }else if(csc$step1$filetype=="csv"){
      pbmc_data=read.table(csc$indata[[i]],sep = csc$step1$csv_sep,header = T,row.names = 1)
    }else{logerror("this is cant open ,please cheack")
      exists()}
    pbmc=CreateSeuratObject(counts = pbmc_data,project = names(csc$indata[i])) #project in orig.ident
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    col2=length(pbmc@active.ident)
 
    mt_vln=VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=-1) 
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_voilin"),sep="/"),mt_vln)
    featureqc1=FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")+ geom_hline(aes(yintercept=csc$step1$percent_mt[2]),colour="#BB0000", linetype="dashed") + ggtitle("QCmt")
    featureqc2=FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_hline(aes(yintercept=csc$step1$nFeature_RNA[1]),colour="#BB0000", linetype="dashed") + geom_hline(aes(yintercept=csc$step1$nFeature_RNA[2]),colour="#BB0000", linetype="dashed") + ggtitle("QCng")
    
    df <- data.frame(nc = pbmc@meta.data$nCount_RNA,pm = pbmc@meta.data$percent.mt)
    p <- ggplot(df,aes(x=nc,y=pm,colour=nc))+geom_point(size=1,shape=16)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+scale_y_continuous(expand=c(0,0))+theme(panel.border = element_blank(), axis.line = element_line(size = 0.7))+scale_color_gradientn(colours = c("#0066CC","#FFCC33","#FF0000"))+labs(x="nCount_RNA",y="percent.mt",colour="nCount_RNA")
    p <- p+theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.x = element_text(size=14),axis.title.y = element_text(size = 14))+geom_hline(aes(yintercept=20),colour="#BB0000", linetype="dashed")
    
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_countVmt"),sep="/"),featureqc1)
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_countVfeature"),sep="/"),featureqc2)
    savepdf(paste(outdir,"step1",names(csc$indata[i]),paste0(names(csc$indata[i]),"_libraryVmt"),sep="/"),p)    

    pbmc<- subset(pbmc,subset = nFeature_RNA > csc$step1$nFeature_RNA[1]& nFeature_RNA < csc$step1$nFeature_RNA[2] & percent.mt >  csc$step1$percent_mt[1] & percent.mt < csc$step1$percent_mt[2])
    
    col3=length(pbmc@active.ident)
    
    listobjall[[i]]=pbmc
    
    cat(pbmc@project.name,
        col2,col3,col2-col3,"\n",
        sep = "\t",
        file = paste(outdir,"step1","summary.xls",sep="/"),fill= F, labels=NULL, append=T)
  }
  loginfo("save %s",paste0(outtemp,"step1.RDS"))
  saveRDS(listobjall,paste0(outtemp,"step1.RDS"))  
  loginfo("step1 end")
}
if (isstep(2)){
  loginfo("step2 start")
  if (!isstep(1)){
    loginfo("input step1 rds")
    listobjall=getrds(1)
  }
  if(csc$step2$normethod=="SCT"){
    filternum=NULL
    for (i in 1:length(listobjall)) {
      listobjall[[i]] <- SCTransform(listobjall[[i]])
      if(csc$step2$kfilter>length(listobjall[[i]]@active.ident)){
        logwarn("%s numcell %s less than %s,plesae check it" ,listobjall[[i]]@project.name,length(listobjall[[i]]@active.ident),csc$step2$kfilter)
        
        next
      }
      filternum=c(filternum,i)
      
    }
    listobj=listobjall[filternum]
    #listname<- unlist(lapply(listobj,levels))
    listname<- unlist(lapply(listobj,function(x){x@project.name}))
    usereflist<-unlist(lapply(csc$step2$reference,function(x){which(listname==x)}))
    if(length(usereflist)==0){
      usereflist=NULL
    }else{    logwarn("userlist in samplelist ")
}
    pbmc_features<-SelectIntegrationFeatures(object.list = listobj, nfeatures = csc$step2$nFeature)
    listobj <- PrepSCTIntegration(object.list = listobj, anchor.features = pbmc_features)
    pbmc_anchors <- FindIntegrationAnchors(object.list = listobj, normalization.method = "SCT", k.filter = csc$step2$kfilter,anchor.features = pbmc_features,reference = usereflist)
    pbmcall <- IntegrateData(anchorset = pbmc_anchors, normalization.method = "SCT")
  }else if(csc$step2$normethod=="vst"){
    filternum=NULL
    for (i in 1:length(listobjall)) {
      logdebug("ok")
      if(csc$step2$kfilter>length(listobjall[[i]]@active.ident)){
        logwarn("%s numcell %s less than %s,plesae check it" ,listobjall[[i]]@project.name,length(listobjall[[i]]@active.ident),csc$step2$kfilter)
        
      }else{filternum=c(filternum,i)}
      listobjall[[i]] <- NormalizeData(listobjall[[i]])
      listobjall[[i]] <- FindVariableFeatures(listobjall[[i]], selection.method = "vst",nfeatures = csc$step2$nFeature)
    }
    listobj=listobjall[filternum]
    pbmc_anchors <- FindIntegrationAnchors(object.list = listobj,normalization.method = "LogNormalize",k.filter = csc$step2$kfilter)
    pbmcall <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)
    pbmcall <- ScaleData(pbmcall, verbose = FALSE)
  }else if(csc$step2$normethod=="none"){
    mergepb=merge(listobjall[[1]],listobjall[-1])
    mergepb=NormalizeData(mergepb,normalization.method = "LogNormalize")
    mergepb=FindVariableFeatures(mergepb,selection.method = "vst",nfeatures = csc$step2$nFeature)
    pbmcall <- ScaleData(mergepb, verbose = FALSE)
  }
  loginfo("save %s",paste0(outtemp,"step2.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step2.RDS"))
  loginfo("step2 end")
}

if (isstep(3)){
  loginfo("step3 start")
  dir.create(paste(outdir,"step3",sep="/"))
  if (!isstep(2)){
    loginfo("input step2 rds")
    pbmcall=getrds(2)
  }
  pbmcall=RunPCA(pbmcall)
  dir.create(paste(outdir,"step3","pca",sep="/"),recursive = T)
  fig1=DimPlot(object = pbmcall, reduction = "pca",group.by = "orig.ident")+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  savepdf(paste(outdir,"step3","pca","pca",sep="/"),fig1,widthle = 10,heightle = 8)
  for(i in 1:csc$step3$heatmapnumber){
    pdf(paste(outdir,"step3","pca",paste0("pcaheatmap",i,".pdf"),sep="/"))
    DimHeatmap(object = pbmcall ,cells=1000, balanced = TRUE ,reduction = "pca",dims=i)
    dev.off()
    png(paste(outdir,"step3","pca",paste0("pcaheatmap",i,".png"),sep="/"))
    DimHeatmap(object = pbmcall ,cells=1000, balanced = TRUE ,reduction = "pca",dims=i)
    dev.off()
  }
  bow_plot=ElbowPlot(pbmcall,csc$step3$elbowdims)#    add parameters: csc$step3$elbowdims by zhangdx
  savepdf(paste(outdir,"step3","pca","bowplot",sep="/"),bow_plot)
  
  cscumap=csc$step3$reduction
  
  if(cscumap=="umap"){
    #pbmcall=RunUMAP(pbmcall,dims = 1:csc$step3$dims,umap.method="umap-learn",metric="correlation")
    pbmcall=RunUMAP(pbmcall,dims = 1:csc$step3$dims,metric="correlation")
    
  }else{
    pbmcall=RunTSNE(pbmcall,dims=1:csc$step3$dims)
  }
  dir.create(paste(outdir,"step3",cscumap,sep="/"))
  umaporig=DimPlot(pbmcall,reduction = cscumap,split.by = "orig.ident",group.by = "orig.ident",ncol= 4)+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  
  savepdf(paste(outdir,"step3",cscumap,"plotby_ident",sep="/"),umaporig,widthle = min(numsap*4,16),heightle = ceiling(numsap/4)*4)
  umaporig2=DimPlot(pbmcall,reduction = cscumap,group.by = "orig.ident")+
    theme(legend.text = element_text(size =8))+
    guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
  
  savepdf(paste(outdir,"step3",cscumap,"plotall_ident",sep="/"),umaporig2,widthle = 10,heightle = 8)
  
  pbmcall <- AddMetaData(pbmcall,pbmcall@reductions$umap@cell.embeddings,col.name = colnames(pbmcall@reductions$umap@cell.embeddings))
  umap <- ggplot(pbmcall@meta.data,aes(x=UMAP_1,y=UMAP_2))+geom_point(size =1,aes(color=log2(nCount_RNA)))+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+theme(panel.border = element_blank(), axis.line = element_line(size = 0.7))+scale_color_gradientn(colours = c("#0066CC","#FFCC33","#FF0000"))
  savepdf(paste(outdir,"step3",cscumap,"plotby_nCount",sep="/"),umap,widthle = 10,heightle = 8)
  
  if(csc$step3$clustercell){
    pbmcall <- FindNeighbors(pbmcall, dims = 1:csc$step3$dims)
    pbmcall <- FindClusters(pbmcall, resolution = csc$step3$resolution,algorithm = csc$step3$algorithm)
    umaporig3=DimPlot(pbmcall,reduction = cscumap,group.by = "seurat_clusters")+
      theme(legend.text = element_text(size =12))+
      guides(color=guide_legend(ncol = 1,override.aes = list(size=3)))
    
    savepdf(paste(outdir,"step3",cscumap,"plotby_cluster",sep="/"),umaporig3,widthle = 10,heightle = 8)
    #pbmcall.markers <- FindAllMarkers(pbmcall,only.pos = TRUE,)
    ##########################################################
    test=pbmcall@meta.data %>% group_by(seurat_clusters,orig.ident) %>% summarise(count=n()) %>% spread(key = seurat_clusters,value = count)
    test2=pbmcall@meta.data %>% group_by(orig.ident) %>% summarise(count=n())
    test3=merge(test2,test)
    test3[is.na(test3)]<-0

    datasum <- apply(test3[,-1],2,sum)
    sumcol <- c("all",datasum)
    finaldata <- rbind(test3,sumcol)
    
    write.table(finaldata,file = paste(outdir,"step3",cscumap,"summary.xls",sep="/"),sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

    total3=test3 %>% gather(cluster,counts,-orig.ident)
    barpli = ggplot(total3,aes(cluster,counts,fill=orig.ident)) + geom_bar(stat = 'identity',position = 'fill') + scale_x_discrete(limits = colnames(test3)[-1])

    savepdf(paste(outdir,"step3",cscumap,"barplot",sep="/"),barpli,widthle = 10, heightle = 8) 
   ###############################################################
    hpca.se <- readRDS("/data/haplox/users/liyq/singlecell/bin/singleRdata/test.rds")
    largdata=readRDS("/data/haplox/users/liyq/singlecell/bin/singleRdata/DatabaseImmuneCellExpressionDataRDS.rds")

    testcd=GetAssayData(pbmcall)
    immannot=SingleR(test = testcd,ref = hpca.se,labels = hpca.se$label.fine)
    foo=immannot$labels
    names(foo)=immannot@rownames
    pbmcall$singleR=foo

    singRtxt= pbmcall@meta.data %>% group_by(seurat_clusters,singleR) %>% summarise(n_cells=n()) %>% arrange(seurat_clusters,desc(n_cells))
    x <- vector()
    for (i in levels(pbmcall$seurat_clusters)) {
    spart = subset(singRtxt,singRtxt$seurat_clusters == i)
    spart = subset(spart,(spart$n_cells/sum(spart$n_cells))>0.02)
    x <- union(x,spart$singleR)
    }  
    x <- sort(x)
    pbmcall1 = pbmcall
    for (j in 1:length(pbmcall$singleR)) {
       if (!pbmcall$singleR[j] %in% x){
       pbmcall1$singleR[j] = 'others'
      }
    }

    cols <- brewer.pal(7,"Spectral")
    pal <- colorRampPalette(cols[c(7,6,2,1)])
    mycolors <- pal(length(x)+1)

    anovplot = DimPlot(pbmcall1,reduction = cscumap,label = F,cols= mycolors,group.by = "singleR")
    levels(anovplot$data$singleR) <- c(x,'others')
    
    anovplot <- anovplot + theme(legend.text = element_text(size =5))+ 
      guides(color=guide_legend(ncol = 1,keyheight = 0.5,override.aes = list(size=1.5))) 
    savepdf(paste(outdir,"step3",cscumap,"celltype",sep="/"),anovplot,widthle = 10, heightle = 8)

    ck <- c(x,'others')
    dzb <- cbind(ck,mycolors)

    write.table(singRtxt,file=paste(outdir,"step3",cscumap,"singleR_celltype.xls",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)

    s1 = pbmcall1@meta.data %>% group_by(seurat_clusters,singleR) %>% summarise(n_cells=n()) %>% arrange(seurat_clusters,desc(n_cells))
    for (i in levels(pbmcall$seurat_clusters)) {

    s10 = subset(s1,s1$seurat_clusters== i)

    s11 = s10 %>% arrange(singleR)
    s11$singleR=factor(s11$singleR,levels=levels(anovplot$data$singleR))
    s11 = s11 %>% arrange(factor(s11$singleR))

    thiscolor <- vector()
    for (j in 1:length(s11$singleR)) {
       for (z in 1:length(dzb[,1])) {
          if (s11$singleR[j]==dzb[z,1]) {
             thiscolor[j] = dzb[z,2]
        }
       }
     } 

    p <- ggplot(data = s11, mapping = aes(x='Content',y = n_cells, fill = singleR)) + geom_bar(stat = 'identity', position = 'stack',width=1,aes(fill = singleR))
    aaa <- round(s11$n_cells/sum(s11$n_cells)*100,2)
    label_value <- ifelse(aaa > 2, paste('(',aaa,'%)',sep=''),"")
    p <- p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +theme(panel.grid.major=element_line(colour=NA),panel.border=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) +labs(fill= paste("Cell type in cluster ",i,sep = ''))
    p <- p+geom_text(stat="identity",aes(x =sum(s11$n_cells)/(sum(s11$n_cells/1.8)),y=s11$n_cells, label = label_value), size=3, position=position_stack(vjust = 0.5))+scale_fill_manual(values = thiscolor)
    
    savepdf(paste(outdir,"step3",cscumap,paste("cluster",i,sep = ''),sep="/"),p,widthle = 10, heightle = 8)
    }  

   
  }
  loginfo("save %s",paste0(outtemp,"step3.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step3.RDS"))
  loginfo("step3 end")
  
}

if (isstep(4)){
  loginfo("step4 start")
  dir.create(paste(outdir,"step4",sep="/"))
  if (!isstep(3)){
    loginfo("input step3 rds")
    pbmcall=getrds(3)
  }
  if(csc$step4$clustermarkers){
    dir.create(paste(outdir,"step4",sep="/"))
    pbmcall.markers <- FindAllMarkers(pbmcall,only.pos = TRUE,test.use = csc$step4$findmarkers_testuse,min.pct = csc$step4$min_pct)
    top10<- subset(pbmcall.markers,avg_logFC > 0.5 & p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, avg_logFC)
    sorttop10 <- top10[order(top10$cluster,-top10$avg_logFC),]
    write.table(sorttop10,file=paste(outdir,"step4",paste0("clusterall_top10genelist.xlsx"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
    write.table(subset(pbmcall.markers,avg_logFC > 0.5 & p_val_adj < 0.05),file=paste(outdir,"step4",paste0("clusterall_adj0.05_logFC0.5genelist.xlsx"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
    
    #pdf(paste(outdir,"step4",paste0("cluster_top10geneheatmap.pdf"),sep="/"))
    #DoHeatmap(pbmcall,features = top10$gene) + scale_y_discrete(breaks=NULL)
    #dev.off()
    plot <- DoHeatmap(object = pbmcall,features = top10$gene) +NoLegend() + scale_y_discrete(breaks=NULL)
    savepdf(paste(outdir,"step4",paste0("cluster_top10geneheatmap"),sep="/"),plot)
    #png(paste(outdir,"step4",paste0("cluster_top10geneheatmap.png"),sep="/"))
    #DoHeatmap(pbmcall,features = top10$gene) + scale_y_discrete(breaks=NULL)
    #dev.off()
    
    pbmcall@active.assay <- "RNA"
    for (i in levels(pbmcall$seurat_clusters)){
      logdebug(i)
      dir.create(paste(outdir,"step4",paste0("cluster",i),sep="/"))
      write.table(subset(pbmcall.markers,cluster==i),file=paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genelist.xlsx"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
      write.table(subset(pbmcall.markers,cluster==i & avg_logFC > 0.5 & p_val_adj < 0.05),file=paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_padj0.05_logFC0.5genelist.xlsx"),sep="/"),sep='\t',row.names = TRUE, col.names = NA,quote = F)
      if(length(subset(top10,cluster==i)$gene)!=0){
        logdebug("abc")
        
        fig4.1=VlnPlot(pbmcall,features=c(subset(top10,cluster==i)$gene),ncol=4,pt.size=-1)
        savepdf(paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genevlnplot"),sep="/"),fig4.1,widthle = 20,heightle = 15)###add l
        fig4.2=FeaturePlot( pbmcall, 
                            features = c(subset(sorttop10,cluster==i)$gene), 
                            cols = c("lightgrey", "blue"),
                            ncol = 4,reduction = csc$step3$reduction)
        savepdf(paste(outdir,"step4",paste0('cluster',i,"/cluster",i,"_genereduction"),sep="/"),fig4.2,widthle = 20,heightle = 15)
      }
    } 
    if(!is.null(csc$step4$custer)){
      logdebug("ok")
      dir.create(paste(outdir,"step4","custer",sep="/"))
      for(i in names(csc$step4$custer)){
        
        igene=eval(parse(text=paste0("csc$step4$custer$",i)))
        dir.create(paste(outdir,"step4","custer",i,sep="/"))
        noexistgene=setdiff(igene,row.names(pbmcall@assays$RNA))
        igene=intersect(igene,row.names(pbmcall@assays$RNA))
        write.table(igene,file = paste(outdir,"step4","custer",i,"use_genelist",sep="/"),sep='\n',row.names = F, col.names = F,quote = F)
        
        logwarn("list of not exist gene: %s" , noexistgene)
        for (m in igene){
          logdebug(m)
          fig4.3=VlnPlot(pbmcall,features = m,pt.size=-1)
          savepdf(paste(outdir,"step4","custer",i,paste0(m,"_vlnplot"),sep="/"),fig4.3)
          fig4.4=FeaturePlot( pbmcall, 
                              features = m, 
                              cols = c("lightgrey", "blue"),
                              reduction = csc$step3$reduction)
          savepdf(paste(outdir,"step4","custer",i,paste0(m,"_reduction"),sep="/"),fig4.4)
        }
        if(length(names(pbmcall@assays))>1){
          pbmcall@active.assay <- names(pbmcall@assays)[2]}
        fig4.5=DotPlot(pbmcall, features = rev(igene), dot.scale = 8) +
          scale_color_gradient(low = "#132B43", high = "#56B1F7") +
          theme(axis.text.x  = element_text(angle = 90))
        savepdf(paste(outdir,"step4","custer",i,"use_genedotplot",sep="/"),fig4.5,widthle = 3+1.5*length(igene),heightle = 4+0.35*length(levels(pbmcall$seurat_clusters)))
      }
      
    }
    
  }
  if(!is.null(csc$step4$difcluster)){
    logdebug("difcluster为寻找不同分组的cluster之间的差异基因")
    dir.create(paste(outdir,"step4","difcluster",sep="/"))
    for (i in names(csc$step4$difcluster)){
      logdebug(i)
      idif=eval(parse(text=paste0("csc$step4$difcluster$",i)))
      dir.create(paste(outdir,"step4","difcluster",i,sep="/"))
      difgenelist=FindMarkers(pbmcall,ident.1 = idif$a,ident.2 = idif$b,group.by = "seurat_clusters",test.use = idif$testuse)
      difgenelist$gene <- rownames(difgenelist)
      dif_filgenelist=subset(difgenelist, avg_logFC > 0.5 & p_val_adj < 0.05)
      write.table(difgenelist,file=paste(outdir,"step4","difcluster",i,"dif_allgenelist.xlsx",sep="/"),sep='\t', row.names = F,col.names = TRUE,quote = F)
      write.table(dif_filgenelist,file=paste(outdir,"step4","difcluster",i,"dif_logFC0.5padj0.05genelist.xlsx",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)
      #展示差异的前十个基因
      arrange(dif_filgenelist,desc(avg_logFC))
      fig4.2=FeaturePlot( pbmcall, 
                          features = head(arrange(dif_filgenelist,desc(avg_logFC)),n = 10)$gene, 
                          cols = c("lightgrey", "blue"),
                          ncol = 4,reduction = csc$step3$reduction)
      savepdf(paste(outdir,"step4","difcluster",i,"dif_genereduction",sep="/"),fig4.2,widthle = 20,heightle = 15)
      if (!is.null(csc$step4$ClusterProfiler) & csc$step4$ClusterProfiler[1]=="true" ){
        dir.create(paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"))
        setwd(paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"))
        cat(csc$step4$ClusterProfiler[-1],
            "-f",paste(outdir,"step4","difcluster",i,"dif_logFC0.5padj0.05genelist.xlsx",sep="/"),
            "-n",paste("\"",csc$title,"\""),
            "-o",paste(outdir,"step4","difcluster",i,"clusterprofiler",sep="/"),
            file = paste(outdir,"step4","difcluster",i,"clusterprofiler","run.sh",sep="/"),
            "&",
            sep=" ", fill= F, labels=NULL, append=F)
        system(paste0("sh ",paste(outdir,"step4","difcluster",i,"clusterprofiler","run.sh",sep="/")))
      }
    }
    
  }
  if(!is.null(csc$step4$difident)){
    logdebug("difident为寻找样品间差异基因")
    dir.create(paste(outdir,"step4","difident",sep="/"))
    for (i in names(csc$step4$difident)){
      logdebug(i)
      idif=eval(parse(text=paste0("csc$step4$difident$",i)))
      dir.create(paste(outdir,"step4","difident",i,sep="/"))
      
      
      difgenelist=FindMarkers(pbmcall,ident.1 = idif$a,ident.2 = idif$b,group.by = "orig.ident",test.use = idif$testuse)
      difgenelist$gene <- rownames(difgenelist)
      dif_filgenelist=subset(difgenelist, avg_logFC > 0.5 & p_val_adj < 0.05)
      write.table(difgenelist,file=paste(outdir,"step4","difident",i,"dif_allgenelist.xlsx",sep="/"),sep='\t', row.names = F,col.names = TRUE,quote = F)
      write.table(dif_filgenelist,file=paste(outdir,"step4","difident",i,"dif_logFC0.5padj0.05genelist.xlsx",sep="/"),sep='\t',row.names = F, col.names = TRUE,quote = F)
      #展示差异的前十个基因
      arrange(dif_filgenelist,desc(avg_logFC))
      fig4.2=FeaturePlot( pbmcall, 
                          features = head(arrange(dif_filgenelist,desc(avg_logFC)),n = 10)$gene, 
                          cols = c("lightgrey", "blue"),
                          ncol = 4,reduction = csc$step3$reduction)
      savepdf(paste(outdir,"step4","difident",i,"dif_genereduction",sep="/"),fig4.2,widthle = 20,heightle = 15)
      if (!is.null(csc$step4$ClusterProfiler) & csc$step4$ClusterProfiler[1]=="true" ){
        dir.create(paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"))
        setwd(paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"))
        cat(csc$step4$ClusterProfiler[-1],
            "-f",paste(outdir,"step4","difident",i,"dif_logFC0.5padj0.05genelist.xlsx",sep="/"),
            "-n",paste("\"",csc$title,"\""),
            "-o",paste(outdir,"step4","difident",i,"clusterprofiler",sep="/"),
            file = paste(outdir,"step4","difident",i,"clusterprofiler","run.sh",sep="/"),
            "&",
            sep=" ", fill= F, labels=NULL, append=F)
        system(paste0("sh ",paste(outdir,"step4","difident",i,"clusterprofiler","run.sh",sep="/")))
      }
    }
    
    
  }
  loginfo("save %s",paste0(outtemp,"step4.RDS"))
  saveRDS(pbmcall,paste0(outtemp,"step4.RDS"))
  loginfo("step4 end")
  
}
if (isstep(5)){
  loginfo("step5 start")
  dir.create(paste(outdir,"step5",sep="/"))
  if (!isstep(4)){
    loginfo("input step3 rds")
    pbmcall=getrds(4)
  }
  library(monocle)
  #从Seurat转换
  datamon=as(as.matrix(pbmcall@assays$RNA@data),'sparseMatrix')
  pd<-new("AnnotatedDataFrame", data = pbmcall@meta.data)
  fd<-new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(datamon), row.names = row.names(datamon)))
  cds <- newCellDataSet(datamon, phenoData = pd, featureData = fd)
  sce = cds
  
  sce <- estimateSizeFactors(sce)
  sce <- estimateDispersions(sce)
  
  #过滤细胞QC
  sce <- detectGenes(sce, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(sce),num_cells_expressed >= 10))
  pData(sce)$Total_mRNAs <- Matrix::colSums(exprs(sce))
  sce <- sce[, pData(sce)$Total_mRNAs < 1e6 ]
  upper_bound <- 10^(mean(log10(pData(sce)$Total_mRNAs)) +
                       2*sd(log10(pData(sce)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(sce)$Total_mRNAs)) -
                       2*sd(log10(pData(sce)$Total_mRNAs)))
  fig5.1=qplot(Total_mRNAs, data = pData(sce), color = orig.ident, geom = "density") +
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound) 
  savepdf(paste(outdir,"step5","QC",sep="/"),fig5.1)
  
  sce <- sce[,pData(sce)$Total_mRNAs > lower_bound &
               pData(sce)$Total_mRNAs < upper_bound]
  sce <- detectGenes(sce, min_expr = 0.1)
  
  
  #选择合适的基因来标记状态
  #diff_test_res <- differentialGeneTest(sce[expressed_genes,],
  #                                     fullModelFormulaStr = "~seurat_clusters")
  #sce_ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  disp_table <- dispersionTable(sce)
  sce_ordering_genes <- subset(disp_table,mean_expression >= csc$step5$meanexpression &dispersion_empirical >= 1 * dispersion_fit)$gene_id
  sce <-
    setOrderingFilter(sce,
                      ordering_genes = sce_ordering_genes)
  fig5.2=plot_ordering_genes(sce)
  dir.create(paste(outdir,"step5","pseudotime",sep="/"))
  savepdf(paste(outdir,"step5","pseudotime","geneplot",sep="/"),fig5.2)
  
  
  sce <- reduceDimension(sce, method = 'DDRTree')
  sce <- orderCells(sce)
  fig5.3=plot_cell_trajectory(sce, color_by = "seurat_clusters")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_byclusters",sep="/"),fig5.3)
  fig5.4=plot_cell_trajectory(sce, color_by = "Pseudotime")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_byPseudotime",sep="/"),fig5.4)
  fig5.5=plot_cell_trajectory(sce, color_by = "seurat_clusters") +     facet_wrap(~seurat_clusters, nrow = ceiling(length(levels(sce$seurat_clusters))/4))
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_splitbyclusters",sep="/"),fig5.5,widthle = 12,heightle = 4+3*ceiling(length(levels(sce$seurat_clusters))/4))
  fig5.6=plot_cell_trajectory(sce, color_by = "orig.ident")+     facet_wrap(~orig.ident, nrow = ceiling(length(union(sce$orig.ident,NULL))/4))
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_splitbyorig.ident",sep="/"),fig5.6,widthle = 12,heightle = 4+3*ceiling(length(union(sce$orig.ident,NULL))/4))
  fig5.7=plot_cell_trajectory(sce, color_by = "State")
  savepdf(paste(outdir,"step5","pseudotime","pseudotime_bystate",sep="/"),fig5.7)
  write.table(select(pData(sce),State,orig.ident,seurat_clusters),file=paste(outdir,"step5","pseudotime","cell_state_clusters.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  
  #差异表达分析
  diff_test_res <- differentialGeneTest(sce[expressed_genes, ],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_genes <- subset(diff_test_res, qval < 0.01)
  sig_genes <- sig_genes[order(sig_genes$qval), ]
  dir.create(paste(outdir,"step5","diffgenes",sep="/"))
  write.table(select(sig_genes,pval,qval,num_cells_expressed),file=paste(outdir,"step5","diffgenes","cell_diffgene.csv",sep="/"),sep='\t',row.names = T, col.names = NA,quote = F)
  
  pdf(paste(outdir,"step5","diffgenes","topgeneheatmap.pdf",sep="/"))
  fig5.8=plot_pseudotime_heatmap(sce[row.names(sig_genes)[1:csc$step5$genenum],],
                                 num_clusters = csc$step5$numclusters,
                                 cores = 1,
                                 show_rownames = T)
  
  dev.off()
  #BEAM分析
  BEAM_res <- BEAM(sce[expressed_genes, ], branch_point = csc$step5$pointid, cores = 16)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  dir.create(paste(outdir,"step5","BEAM",sep="/"))
  write.table(BEAM_res,file=paste(outdir,"step5","BEAM","diffgenelist.csv",sep="/"),sep='\t',row.names = F, col.names = T,quote = F)
  pdf(paste(outdir,"step5","BEAM","topgeneheatmap.pdf",sep="/"))
  
  plot_genes_branched_heatmap(sce[row.names(subset(BEAM_res,
                                                   qval < 1e-4))[1:csc$step5$BEAMgn],],
                              branch_point = csc$step5$pointid,
                              num_clusters = csc$step5$BEAMnumclusters,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T)
  
  
  dev.off()
  fig5.9=plot_genes_branched_pseudotime(sce[rownames(subset(BEAM_res, qval < 1e-8))[1:9], ], 
                                        branch_point = csc$step5$pointid, color_by = "State", 
                                        ncol = 3)
  savepdf(paste(outdir,"step5","BEAM","top10gene",sep="/"),fig5.9)
  
  if(!is.null(csc$step5$BEAMgenelist)){
    
    fig5.10=plot_genes_branched_pseudotime(sce[intersect( rownames(sce),csc$step5$BEAMgenelist),],branch_point = csc$step5$pointid, color_by = "State", 
                                           ncol = 3)
    savepdf(paste(outdir,"step5","BEAM","custergene",sep="/"),fig5.10)}
  
  loginfo("save %s",paste0(outtemp,"step5.RDS"))
  saveRDS(sce,paste0(outtemp,"step5.RDS"))
  loginfo("step5 end")
  
  
  
}



isstep(1)
if(csc$run$step1==F){
  loginfo("ok")
  
  csc=parseTOML(opt$infile)
  pbmcall=RunPCA(pbmcall1)
  pbmcall1=RunUMAP(pbmcall1,dims = 1:30,umap.method="umap-learn",metric="correlation")
  DimPlot(pbmcall1,split.by = "orig.ident")
  VlnPlot(pbmcall,features=c(subset(top10,cluster==4)$gene),ncol=4,pt.size=-1)
  FeaturePlot( pbmcall, 
               features = c(subset(sorttop10,cluster==5)$gene), 
               cols = c("lightgrey", "blue"), reduction = "umap")
  DoHeatmap(pbmcall,features = top10$gene)
  write.table(nor_gsva_matrix1,file=paste(outdir,"step4",'nor_gsva_matrix.xls',sep="/"),sep='\t',row.names = TRUE, col.names = TRUE,quote = F)
  top10 <- subset(pbmcall.markers,avg_logFC > 1 & p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, avg_logFC)
}
