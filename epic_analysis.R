suppressPackageStartupMessages({
    library("R.utils")
    library("missMethyl")
    library("limma")
    library("minfi")
    library("IlluminaHumanMethylationEPICmanifest")
    library("RColorBrewer")
    library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    library("eulerr")
    library("plyr")
    library("gplots")
    library("reshape2")
    library("beeswarm")
})

# probably need some annotations
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# import dataset
baseDir <- "ILMLEPIC-15799"
targets <- read.metharray.sheet(baseDir)
head(targets)
rgSet <- read.metharray.exp(targets = targets)
rgSet

# normalisation
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=FALSE)
par(mfrow=c(1,2), cex=.7)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")

# filtering
# include sex chromosomes
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]
# exclude sex chromosomes
keep <- !(featureNames(mSetSw) %in% annEPIC$Name[annEPIC$chr %in% 
                                                   c("chrX","chrY")])
mSetFlt <- mSetSw[keep,]

# normalisation include sex chromosomes
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)
dim(Mval)

# normalisation exclude sex chromosomes
meth <- getMeth(mSetFlt)
unmeth <- getUnmeth(mSetFlt)
Mval_flt <- log2((meth + 100)/(unmeth + 100))
beta_flt <- getBeta(mSetFlt)
dim(Mval_flt)

# MDS
par(mfrow=c(1,2), cex=.7)
colour_palette=brewer.pal(n = length(levels(targets$Sample_Group)), name = "Paired")
colors <- colour_palette[as.integer(factor(targets$Sample_Group))]
plotMDS(Mval, labels=mylabels,col=colors,main="sex chromosomes included")
legend("bottomright",legend=levels(factor(targets$Sample_Group)),pch=16,cex=1,col=colour_palette)
plotMDS(Mval_flt, labels=mylabels,col=colors,main="sex chromosomes excluded")
