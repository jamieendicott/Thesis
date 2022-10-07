library(transcripTools)
library(dplyr)
library(ggplot2)

setwd('/secondary/projects/laird/jamie/0922.EPIC/LAIP_20220826_EPIC/')
samples<-read.csv('0822.samples.csv')
s<-samples[,c(6,12,14,16)]
s<-subset(s,s$experiment=='HGPS' | s$experiment=='BloomSyndrome' | s$experiment=='Werner1' | s$experiment=='Werner2')
betas<-read.csv('JP.betas.0822.csv',row.names=1,check.names=F)
b<-betas[,c(match(s$EPIC_ID,colnames(betas)))]

#pull in PCGT annotation
pcgt<-read.table('20180904_EPIC.probe.cg.REMC.FIB.E055.H3K27me3.rda', header=TRUE, quote="\"")

pb<-subset(b,rownames(b)%in%pcgt$RDX2)
dim(pb)
#[1] 226833     25
s$pcgt<-apply(pb,2,median,na.rm=t)
#weirdly these lose methylation... so many

#pull mostvar 5k, pheatmap
library(pheatmap)
library(viridisLite)
library(transcripTools)
pb<-na.omit(pb)
dat<-mostVar(pb3,5000)
anno<-s[,c(1,4)]
anno$experiment<-as.factor(anno$experiment)
rownames(anno)<-s$EPIC_ID

#b2<-b[order(b[,1],decreasing = TRUE),]

p<-pheatmap(dat,cluster_cols = T, cluster_rows =T,
         show_rownames = F, show_colnames = F,
         color=turbo(100), annotation_col=anno2
)

pdf('pcgt.ctl.hm.pdf')
p
dev.off()
