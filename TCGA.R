read.table("/haplox/haprs/wuyuling/project/HPLS2019050602J_liver_20190710/liyq/TCGA/database/TCGA-LIHC.survival.tsv")
library(VennDiagram)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggThemeAssist)
setwd("/haplox/haprs/liyq/other/TCGGA/TCGA")
fpkmdata=read.table("database/name_symbol_gene_fpkm.tsv",header = T,stringsAsFactors = F)
survialdata=read.table("database/name_survival.tsv",header = T,stringsAsFactors = F)
vapheno=read.table("tempdata/va_pheno",header = T, stringsAsFactors = F)
nor_desk=read.table("TIL/IPS/desk.xls",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
new_phone_all=read.table("/haplox/haprs/wuyuling/project/HPLS2019050602J_liver_20190710/liyq/TCGA/TIL/all_group.txt",sep = "\t",header = T,row.names = 1)
outpath="/haplox/haprs/liyq/other/TCGGA/TCGA/TIL/IPS/hubo_singcle"
outfile=paste0(outpath,"/","venn.tiff")
rownames(survialdata)=survialdata[,1]
rownames(vapheno)=vapheno[,1]
xxx=Reduce(intersect,list(vapheno[,1],survialdata[,1],colnames(fpkmdata)))

c=list(pheno=vapheno[,1],survia=survialdata[,1],fpkm=colnames(fpkmdata))
venn.diagram(c,  alpha=c(0.5,0.5,0.5),
             col = "black",
             lty = "dotted",
             lwd = 3,
             cex = 1.5,
             fontfamily = "serif",
             fontface = "bold",
             cat.cex = 1.5,
   #          cat.pos =c(-10,10),
             cat.fontfamily = "serif",
             filename=outfile)
outgene=get.venn.partitions(c,keep.elements = F)

uniqsam_name=str_subset(xxx,"_11.$",negate = T)
newfpkm=fpkmdata[,c(colnames(fpkmdata)[1],uniqsam_name)]
newsurvial=survialdata[uniqsam_name,]
newvapheno=vapheno[uniqsam_name,]
dnew=list(pheno=newvapheno[,1],survia=newsurvial[,1],fpkm=colnames(newfpkm))
venn.diagram(dnew,  alpha=c(0.5,0.5,0.5),
             col = "black",
             lty = "dotted",
             lwd = 3,
             cex = 1.5,
             fontfamily = "serif",
             fontface = "bold",
             cat.cex = 1.5,
             #          cat.pos =c(-10,10),
             cat.fontfamily = "serif",
             filename=outfile)

t_nor_desk=as.data.frame(t(nor_desk))
ggplot(t_nor_desk, aes(x="0",y = Act.CD8)) + 
  geom_violin() +
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(fill = NA, width = 0.1, colour = "grey")+
  
  theme_bw()

new_nor_desk=nor_desk[,newsurvial[,1]]
new_ph=pheno[uniqsam_name,,drop= F]
new_ph
testhph=head(new_ph)
new_ph$name=rownames(new_ph)

testhsu=head(newsurvial)
new_phone=merge(new_ph,newsurvial,by.x="name",by.y="sample")
new_phone=merge(new_phone,newvapheno,by.x="name",by.y="submitter_id.samples")
#testhsu[testhph]
phone_all=new_phone[,-c(1,2,4)]

rownames(phone_all)=new_phone[,1]
phone_all[,1]=as.character(phone_all[,1])
listcolors=getlistcolors(phone_all)

outheat=pheatmap(new_nor_desk,show_colnames = T,cluster_rows = F,cluster_cols = clusterc,
         annotation_col = new_phone_all,
         annotation_colors = listcolors$ph,
         #cellwidth=10,cellheight=10,
         filename = paste(outpath,'heatmap_sort.pdf',sep="/"),height = heightdd,width = widthdd)
xxxgroup=cutree(outheat$tree_col,k = 2)
pheno=as.data.frame(xxxgroup)
pheno[,1]=as.character(pheno[,1])
new_phone_all=merge(phone_all,pheno,by = "row.names")
rownames(new_phone_all)=new_phone_all$Row.names
new_phone_all$Row.names=NULL
listcolors=getlistcolors(new_phone_all)
new_phone_all$OS_media=new_phone_all$OS.time>median(new_phone_all$OS.time)

write.table(new_phone_all,file="/haplox/haprs/wuyuling/project/HPLS2019050602J_liver_20190710/liyq/TCGA/TIL/all_group.txt",
            quote =F, row.names = T,col.names = T,sep = "\t")
t_new_nor=t(new_nor_desk)
all_new_nor_h=merge(t_new_nor,new_phone_all,by="row.names")
rownames(all_new_nor_h)=all_new_nor_h$Row.names
all_new_nor_h$Row.names=NULL

ggplot(all_new_nor_h, aes_string(x="vascular_tumor_cell_type" ,y = "MDSC",colour="OS")) + 
  geom_violin() +
  #  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  #geom_jitter(aes(x=vascular_tumor_cell_type,y=MDSC,colour=OS),shape=16, position = position_jitter(width = .05), alpha = 0.5)+
  geom_point(position = position_jitterdodge(dodge.width = 0.9))+
  theme_bw()

all_new_nor=all_new_nor_h
colnames(all_new_nor)
colnames(all_new_nor)=str_replace(colnames(all_new_nor)," ","")

colnames(all_new_nor)=str_replace_all(colnames(all_new_nor)," ","_")


testall=gather(all_new_nor,"celltype","n",1:28)

testall$OS_media=factor(testall$OS_media,levels = c(TRUE,FALSE),ordered = T)
a=ggplot(testall, aes(x=OS_media ,y = n,colour=OS_media)) + 
  geom_violin() +
  geom_boxplot(width=.5)+
  #  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  #geom_jitter(aes(x=vascular_tumor_cell_type,y=MDSC,colour=OS),shape=16, position = position_jitter(width = .05), alpha = 0.5)+
  geom_point(position = position_jitterdodge(jitter.width=1,dodge.width = 0.5),size=0.1,alpha = 0.5)+
  scale_color_manual(values =  c("#5da0f8","#fc7f21"),labels = c("Good Prognosis","Poor Prognosis"))+
  scale_x_discrete(labels=c("Good Prognosis","Poor Prognosis"))+
  ylim(NA,1.1)+
  stat_compare_means(aes(x=OS_media,y=n),method = "wilcox.test",label = "p.signif",label.x.npc="centre",label.y.npc = 0.95)+
  #scale_x_discrete(limits= c("TRUE","FALSE"))+
  theme_bw()+

  facet_wrap(vascular_tumor_cell_type ~ celltype,nrow = 8)
  #facet_wrap(~ vascular_tumor_cell_type, nrow = 6)
savepdf(paste0(outpath,"/plotby_split_celltype_OS"),a,widthle = 14,heightle = 22 )
ggThemeAssistGadget(a)

ex5$memorypar2 <- as.numeric(ex5$tailindex) + 
  3 * (as.numeric(as.character(ex5$memorypar)) - 0.2) 
dodge <- position_dodge(width = 0.9)

p <- ggplot(ex5,aes(x=tailindex , y=hillest)) +
  scale_x_discrete() +
  geom_boxplot(aes(colour = memorypar), outlier.colour = NA, position = dodge) +

  geom_jitter(aes(colour = memorypar, x = memorypar2), 
              position = position_jitter(width = .05), alpha = 0.5) 
  facet_wrap(~ process, nrow = 2)
p


