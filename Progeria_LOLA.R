library(transcripTools)
library(dplyr)
library(ggplot2)
library(LOLA)
library(GenomicRanges)
library(simpleCache)

setwd('/secondary/projects/laird/jamie/0922.EPIC/LAIP_20220826_EPIC/')
betas<-read.csv('JP.betas.0822.csv',row.names=1,check.names=F)
samples<-read.csv('0822.samples.csv')
s<-samples[,c(6,12,14,16)]
s<-subset(s,s$experiment=='HGPS')
b<-betas[,c(match(s$EPIC_ID,colnames(betas)))]

b$delta<-b[,ncol(b)]-b[,1]
hi.me<-subset(b,b[,1]>0.7) 
hi.me<-subset(hi.me,abs(hi.me$delta)<0.1) 
dim(hi.me) #233194
lo.me<-subset(b,b[,1]<0.3)
lo.me<-subset(lo.me,abs(lo.me$delta)<0.1)
dim(lo.me) #209554

regionDB<-loadRegionDB('regions/LOLACore/hg19/') #hg19 has CTCF annotated
man<-read.csv('/secondary/projects/laird/jamie/manifests/EPIC.hg19.manifest.simple.csv',header=T,row.names=1)
U<-man[,1:3]
#convert to granges object
userUniverse<-GRanges(seqnames = U$CpG_chrm,
                      ranges = IRanges(start = U$CpG_beg,
                                       end = U$CpG_end))
loss<-subset(b,b$delta<=-0.5)
group<-lo.me #etc
input<-subset(man,man$probeID%in%rownames(loss)) 
input<-input[,1:3]
userSets<-GRanges(seqnames = input$CpG_chrm,
                  ranges = IRanges(start = input$CpG_beg,
                                   end = input$CpG_beg))
locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)
head(locResults)

writeCombinedEnrichment(locResults, outFolder="loss.lolaRes")

#depletion of DNAm at CTCF sites?
#can plot longitudinally as well, take median. see if wild drop happens?
#first look at delta

CTCF<-read.table('/secondary/projects/laird/jamie/manifests/EPIC.hg19.CTCFbinding.tsv.gz',header=T)
b$CTCF<-CTCF[c(match(rownames(b),CTCF$probeID)),6]
summary(b$CTCF)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      0    1094    5275   12857   15275 1346332   24089
 t.test(loss$CTCF,b$CTCF)

	Welch Two Sample t-test

data:  loss$CTCF and b$CTCF
t = 30.235, df = 6642.7, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 12539.07 14277.75
sample estimates:
mean of x mean of y
 26265.37  12856.95
#CTCF sites don't look to be involved in methylation drop. if anything probably protected?
