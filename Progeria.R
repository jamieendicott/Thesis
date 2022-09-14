#Looking at/reanalyzing preexisting progeria DNAm datasets
#Blood dataset GSE182991 (obviously want to deconvolute samples...)
#Fibroblasts GSE149960
library(GEOquery)
library(GenomicRanges)
library(FlowSorted.Blood.EPIC)
library(rtracklayer)
library(glmnet)
library(coefplot)
library(reshape2)
library(ggplot2)
library(dplyr)
library(sesame)
library(purrr)

setwd('/secondary/projects/laird/jamie/deconvolution/')

geo<-"GSE182991"
getGEOSuppFiles(GEO=geo, makeDirectory = TRUE, baseDir = getwd(),
  fetch_files = TRUE) 
list.files(".")

setwd(paste0('/secondary/projects/laird/jamie/deconvolution/',geo))

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO(geo)
p<-pData(g[[1]])
#trim as needed 
p<-p[,c(1,2,34:36)]
#NOTE: FOR BLOOD ONLY
sset<-SigSetListFromPath(path = ".", parallel = FALSE, recursive = TRUE)
RGset<-SigSetsToRGChannelSet(sset)
data("IDOLOptimizedCpGs") 
head(IDOLOptimizedCpGs)

propEPIC<-estimateCellCounts2(RGset, compositeCellType = "Blood", 
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL", 
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", 
                                "Mono", "Neu"))
percEPIC<-as.data.frame(round(propEPIC$prop*100,1))
percEPIC  
#total myeloid and lymphoid lineages
percEPIC$L<-percEPIC$CD8T+percEPIC$CD4T+percEPIC$NK+percEPIC$Bcell
percEPIC$M<-percEPIC$Mono+percEPIC$Neu
rows<-(as.character(map(strsplit(rownames(percEPIC), split = "_"), 1)))
rownames(percEPIC)<-rows
#add results to p
p<-cbind(p,percEPIC)
write.csv(p,paste0(geo,'.p.csv'))

idat_dir = ('.')
betas = do.call(cbind, lapply(searchIDATprefixes(idat_dir), function(pfx) {
  getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(pfx),pval.threshold = 0.1)))) 
}))
dim(betas)
cols<-(as.character(map(strsplit(colnames(betas), split = "_"), 1)))
colnames(betas)<-cols
#measure median global methylation, median PMD methylation, median PMD solo-WCGW methylation
p$medglobal<-apply(betas,2,median,na.rm=T)
setwd('/secondary/projects/laird/jamie/deconvolution/')
EPIC.comPMD<-read.delim('../manifests/EPIC.comPMD.probes.tsv',header=F)
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)
dim(b)
p$medsoloWCGW<-apply(b,2,median,na.rm=T)
#measure clocks
#Horvath age estimates (sesame funciton)
p$Horv.age<-apply(betas,2,predictAgeHorvath353)
p$Skinblood.age<-apply(betas,2,predictAgeSkinBlood)

#RT
RT<-read.csv('RepliTali_coefs.csv')
x<-betas[c(match(RT$Coefficient[-1],rownames(betas))),c(match(rownames(p),colnames(betas)))]
X<-apply(x,2,function(x) RT$Value[-1]*x)
est.RT<-apply(X,2,function(x) sum(x, na.rm=T) +RT$Value[1])
p$RT<-est.RT
