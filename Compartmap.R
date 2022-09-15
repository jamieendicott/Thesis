library(compartmap)
library(GenomicRanges)
library(Homo.sapiens)
library(sesame)

sset<-SigSetListFromPath(path = ".", parallel = FALSE, recursive = TRUE)
RGset<-SigSetsToRGChannelSet(sset)
b<-sesamize(RGSet)
atac_compartments <- getCompartments(b, type = "atac", parallel = FALSE, chrs = "chr14")
