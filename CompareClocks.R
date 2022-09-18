#for cultured cells etc
library(GEOquery)
library(GenomicRanges)
library(FlowSorted.Blood.EPIC)
library(rtracklayer)
library(glmnet)
library(coefplot)
library(reshape2)
library(ggplot2)
library(sesame)
library(purrr)
#command line edit

#external sets
setwd('/secondary/projects/laird/jamie/GEO/')
geo<-"GSE131280"
getGEOSuppFiles(GEO=geo, makeDirectory = TRUE, baseDir = getwd(),
  fetch_files = FALSE)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO(geo)
p<-pData(g[[1]])
#trim p as needed
p<-p[,c(1,2,34:36)]

setwd(paste0('/secondary/projects/laird/jamie/GEO/',geo))
write.csv(p,paste0(geo,'.p.csv'))

#sesame beta processing
idat_dir = ('.')
betas = do.call(cbind, lapply(searchIDATprefixes(idat_dir), function(pfx) {
  getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(pfx),pval.threshold = 0.1)))) 
}))
dim(betas)

cols<-(as.character(map(strsplit(colnames(betas), split = "_"), 1)))
colnames(betas)<-cols
#Horvath age estimates (sesame funciton)
p$Horv.age<-apply(betas,2,predictAgeHorvath353)
p$Skinblood.age<-apply(betas,2,predictAgeSkinBlood)
#measure median global methylation, median PMD methylation, median PMD solo-WCGW methylation
p$medglobal<-apply(betas,2,median,na.rm=T)
setwd('/secondary/projects/laird/jamie/deconvolution/')
EPIC.comPMD<-read.delim('../manifests/EPIC.comPMD.probes.tsv',header=F)
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)
dim(b)
p$medsoloWCGW<-apply(b,2,median,na.rm=T)

#Apply RepliTali
download.file('https://raw.githubusercontent.com/jamieendicott/Nature_Comm_2022/main/RepliTali/RepliTali_coefs.csv','./RepliTali_coefs.csv')
RT<-read.csv('RepliTali_coefs.csv')
x<-betas[c(match(RT$Coefficient[-1],rownames(betas))),c(match(rownames(p),colnames(betas)))]
X<-apply(x,2,function(x) RT$Value[-1]*x)
est.RT<-apply(X,2,function(x) sum(x, na.rm=T) +RT$Value[1])
p$RT<-est.RT


#epitoc code
epiTOC2 <- function(data.m,ages.v=NULL){
    load("dataETOC2.Rd"); ## this loads the CpG information
    cpgETOC.v <- dataETOC2.l[[2]];
    estETOC2.m <- dataETOC2.l[[1]];
    soloCpG.v <- dataETOC2.l[[3]];
    ### do epiTOC
    common.v <- intersect(rownames(data.m),cpgETOC.v);
    print(paste("Number of represented epiTOC CpGs (max=385)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
    ### do epiTOC2
    map.idx <- match(rownames(estETOC2.m),rownames(data.m));
    rep.idx <- which(is.na(map.idx)==FALSE);
    print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
    tmp.m <- data.m[map.idx[rep.idx],];
    TNSC.v <- 2*colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
    TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
    ### do HypoClock
    common.v <- intersect(rownames(data.m),soloCpG.v);
    print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    hypoSC.v <- colMeans(data.m[map.idx,],na.rm=TRUE);

    estIR.v <- NULL; estIR2.v <- NULL;
    estIR <- NULL;  estIR2 <- NULL;
    if(!is.null(ages.v)){
      estIR.v <- TNSC.v/ages.v;
      estIR <- median(estIR.v,na.rm=TRUE);
      estIR2.v <- TNSC2.v/ages.v;
      estIR2 <- median(estIR2.v,na.rm=TRUE);
    }

    
    return(list(tnsc=TNSC.v,tnsc2=TNSC2.v,irS=estIR.v,irS2=estIR2.v,irT=estIR,irT2=estIR2,pcgtAge=pcgtAge.v,hypoSC=hypoSC.v));
}
setwd('/secondary/projects/laird/jamie/all.EPIC/processed/clock/')
mat<-as.matrix(na.omit(betas))
epi<-epiTOC2(mat)

p$pcgtAge<-epi$pcgtAge
p$hypoSC<-epi$hypoSC
p$epiTOC2.tnsc<-epi$tnsc
p$epiTOC2.tnsc2<-epi$tnsc2

#epiCMIT 
#Bcell specific, prognostic for lymphoid malignancies, unsure if crosses over ?
#code here https://duran-ferrerm.github.io/Pan-B-cell-methylome/Estimate.epiCMIT.html
betas<-as.data.frame(betas)
library(GenomicRanges)
library(data.table)
download.file("https://github.com/Duran-FerrerM/Pan-B-cell-methylome/raw/master/data/Estimate.epiCMIT.RData", destfile = "Estimate.epiCMIT.RData", method="libcurl")
load("Estimate.epiCMIT.RData")

DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = na.omit(betas),
                                DNAm.genome.assembly = "hg19",
                                map.DNAm.to = "Illumina.450K.epiCMIT",
                                min.epiCMIT.CpGs = 800 
)

#1264 from 1348 epiCMIT-CpGs are present in your matrix! ...
epiCMIT.Illumina <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                            return.epiCMIT.annot = FALSE,
                            export.results = FALSE,
                            export.results.dir = ".",
                            export.results.name = "epiCMIT"
)
epiCMIT.Illumina.results <- cbind(epiCMIT.Illumina$epiCMIT.scores,
                                  epiCMIT.CpGs=epiCMIT.Illumina$epiCMIT.run.info$epiCMIT.CpGs,
                                  epiCMIT.hyper.CpGs=epiCMIT.Illumina$epiCMIT.run.info$epiCMIT.hyper.CpGs,
                                  epiCMIT.hypo.CpGs=epiCMIT.Illumina$epiCMIT.run.info$epiCMIT.hypo.CpGs
)

p$epiCMIT<-epiCMIT.Illumina.results$epiCMIT

#MiAge
#scripts available here http://www.columbia.edu/~sw2206/softwares.htm
setwd('../mitotic_age_R_code_new/')
source("function_library.r") ### library of all functions used in the calculation of mitotic ages                                                                                 
clocksitesID=as.matrix(read.csv("Additional_File1.csv",header=T))[,1]   ### epigenetic clock CpG sites
beta=  betas[ match(clocksitesID,rownames(betas)),]##select clock CpG sites 
#268 cpgs
load("site_specific_parameters.Rdata") #### site-specific parameters
b=methyl.age[[1]];c=methyl.age[[2]];d=methyl.age[[3]]
n=mitotic.age(beta,b,c,d) ### estimated mitotic age
names(n)=colnames(beta)

p$MiAge<-n

setwd(paste0('/secondary/projects/laird/jamie/GEO/',geo))
write.csv(p,paste0(geo,'res.all.clocks.csv'))
         
         
