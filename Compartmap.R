library(compartmap)
library(GenomicRanges)
library(Homo.sapiens)
library(sesame)

sset<-SigSetListFromPath(path = ".", parallel = FALSE, recursive = TRUE)
RGset<-SigSetsToRGChannelSet(sset)
b<-sesamize(RGSet)
array_compartments <- getCompartments(ss2, type = "array", parallel = FALSE, chrs = "chr14")

