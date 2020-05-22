

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

# need some annotations
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
myann <- data.frame(annEPIC[,c("UCSC_RefGene_Name","Regulatory_Feature_Group")])
promoters <- grep("Prom",myann$Regulatory_Feature_Group)

# import dataset
baseDir <- "ILMLEPIC-15799"
targets <- read.metharray.sheet(baseDir)
head(targets)
rgSet <- read.metharray.exp(targets = targets)
rgSet

# define sample sheet
mygroups <- read.csv("Groups.csv")
samplesheet <- mygroups[,9:11]
colnames(samplesheet)[1] <- "Sample_Name"

# define control and trt groups for contrast 1
ctrl1_grp <- as.character(mygroups[,1])
ctrl1_grp <- ctrl1_grp[which(ctrl1_grp != "")]
samplesheet$ctrl1 <- as.numeric(samplesheet$Sample_Name %in% ctrl1_grp)
trt1_grp <- as.character(mygroups[,2])
trt1_grp <- trt1_grp[which(trt1_grp != "")]
samplesheet$trt1 <- as.numeric(samplesheet$Sample_Name %in% trt1_grp)

# define control and trt groups for contrast 2
ctrl2_grp <- as.character(mygroups[,5])
ctrl2_grp <- ctrl2_grp[which(ctrl2_grp != "")]
samplesheet$ctrl2 <- as.numeric(samplesheet$Sample_Name %in% ctrl2_grp)
trt2_grp <- as.character(mygroups[,6])
trt2_grp <- trt2_grp[which(trt2_grp != "")]
samplesheet$trt2 <- as.numeric(samplesheet$Sample_Name %in% trt2_grp)

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
keep <- !(featureNames(mSetSw) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetFlt <- mSetSw[keep,]

# normalisation include sex chromosomes
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
colnames(Mval) <- samplesheet$Sample_Name
beta <- getBeta(mSetSw)
dim(Mval)

# normalisation exclude sex chromosomes
meth <- getMeth(mSetFlt)
unmeth <- getUnmeth(mSetFlt)
Mval_flt <- log2((meth + 100)/(unmeth + 100))
colnames(Mval_flt) <- samplesheet$Sample_Name
beta_flt <- getBeta(mSetFlt)
dim(Mval_flt)

# MDS
par(mfrow=c(1,2), cex=.7)
colour_palette=brewer.pal(n = length(levels(targets$Sample_Group)), name = "Paired")
colors <- colour_palette[as.integer(factor(targets$Sample_Group))]
plotMDS(Mval, labels=mylabels,col=colors,main="sex chromosomes included")
legend("bottomright",legend=levels(factor(targets$Sample_Group)),pch=16,cex=1,col=colour_palette)
plotMDS(Mval_flt, labels=targets$Sample_Name,col=colors,main="sex chromosomes excluded")

# DM analysis contrast 1 include sex chromosomes

ss <- subset(samplesheet,ctrl1==1 | trt1==1)
groups <- factor(ss$trt1)
sex <- factor(ss$Gender,levels=c("M","F"))
age <- ss$Age
design <- model.matrix(~ age + sex + groups)
mxs <- Mval[,which( colnames(Mval) %in% ss$Sample_Name )]
fit.reduced <- lmFit(mxs,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
dm1in <- topTable(fit.reduced,coef=4, number = Inf)
dm1ina <- merge(myann,dm1in,by=0)
dm1ina <- dm1ina[order(dm1ina$P.Value),]
head(dm1ina)

# DM analysis contrast 1 excl sex chromosomes
ss <- subset(samplesheet,ctrl1==1 | trt1==1)
groups <- factor(ss$trt1)
sex <- factor(ss$Gender,levels=c("M","F"))
age <- ss$Age
design <- model.matrix(~ age + sex + groups)
mxs <- Mval_flt[,which( colnames(Mval_flt) %in% ss$Sample_Name )]
fit.reduced <- lmFit(mxs,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
dm1ex <- topTable(fit.reduced,coef=4, number = Inf)
dm1exa <- merge(myann,dm1ex,by=0)
dm1exa <- dm1exa[order(dm1exa$P.Value),]
head(dm1exa)

# DM analysis contrast 2 include sex chromosomes
ss <- subset(samplesheet,ctrl2==1 | trt2==1)
groups <- factor(ss$trt2)
sex <- factor(ss$Gender,levels=c("M","F"))
age <- ss$Age
design <- model.matrix(~ age + sex + groups)
mxs <- Mval[,which( colnames(Mval) %in% ss$Sample_Name )]
fit.reduced <- lmFit(mxs,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
dm2in <- topTable(fit.reduced,coef=4, number = Inf)
dm2ina <- merge(myann,dm2in,by=0)
dm2ina <- dm2ina[order(dm2ina$P.Value),]
head(dm2ina)

# DM analysis contrast 2 excl sex chromosomes
ss <- subset(samplesheet,ctrl2==1 | trt2==1)
groups <- factor(ss$trt2)
sex <- factor(ss$Gender,levels=c("M","F"))
age <- ss$Age
design <- model.matrix(~ age + sex + groups)
mxs <- Mval_flt[,which( colnames(Mval_flt) %in% ss$Sample_Name )]
fit.reduced <- lmFit(mxs,design)
fit.reduced <- eBayes(fit.reduced)
summary(decideTests(fit.reduced))
dm2ex <- topTable(fit.reduced,coef=4, number = Inf)
dm2exa <- merge(myann,dm2ex,by=0)
dm2exa <- dm2exa[order(dm2exa$P.Value),]
head(dm2exa)
