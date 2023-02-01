## Universial script for running Essentials analyses
## Arg[1] = normalization method: "TMM", "TMMwsp", "RLE", "upperquantile", "none"
## Arg[2] = p-value adjustment: 'BH', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BY', 'none
## Arg[3] = dispersion: "tagwise", "common"
## Arg[4] = smoothing / prior n: integer
## Arg[5] = analysis type: "qCML", "CR"
## Arg[6] = data type: "gene", "essentialgenes", "ta"

Args <- commandArgs(TRUE)
options(error=expression(NULL))
library(edgeR)

norm.type <- Args[1]
adj.type <- Args[2]
disp.type <- Args[3]
smooth.val <- Args[4]
stat.type <- Args[5]
data.type <- Args[6]

targets <- read.delim(file = paste(data.type, "_targets.txt", sep=""), stringsAsFactors = FALSE) 
d <- readDGE(targets, columns=c(2,1))
d$counts <- round(d$counts)
if (data.type == "gene") d <- d[rowSums(d$counts) > 10*nrow(targets)-1, ]
if (data.type == "essentialgenes") d <- d[rowSums(d$counts) > 10, ]
counts.old <- d$counts
d$counts <- d$counts+1
## !!!! The following step is not originally performed by Essentials (not sure why).
## It updates the library size after the rounding and filtering step.
d$samples$lib.size <- colSums(d$counts)

if (data.type == "gene") title.text = "genes"
if (data.type == "essentialgenes") title.text = "measured and expected counts"
if (data.type == "ta") title.text = "ta sites"
if (nrow(targets) >2){
    p3 <- prcomp(d$counts)
    png(paste(data.type, '_PCA_non_normalized.png', sep=""), width=1280, height=1024)
    biplot(p3, cex=1:2, col=c(40,1), main=paste('PCA Plot non normalized on',title.text, sep=" "))
    dev.off()
    png(paste(data.type,'_MDS_non_normalized.png', sep=""), width=960, height=768)
    plotMDS(d, col=as.numeric(d$samples$group), main=paste("MDS Plot non normalized on", title.text, sep=" "))
    legend("bottomleft", as.character(unique(d$samples$group)), col=1:nlevels(d$samples$group), pch=20)
    dev.off()
} else {
    png(paste(data.type, '_PCA_non_normalized.png', sep=""), width=1280, height=1024)
    plot(log(d$counts))
    title("three samples needed for PCA. x-y plot generated")
    dev.off()
}
                
if (data.type == "essentialgenes"){
    d <- calcNormFactors(d, refColumn=nrow(targets)-1, method=norm.type)
} else {
    d <- calcNormFactors(d, method=norm.type) ## for both gene and ta
}
d <- estimateCommonDisp(d) 

if (nrow(targets) >2){
    p3 <- prcomp(d$pseudo.counts)
    png(paste(data.type,'_PCA_normalized.png', sep=""), width=1280, height=1024)
    biplot(p3, cex=1:2, col=c(40,1), main=paste('PCA Plot normalized on', title.text, sep=" "))
    dev.off()
    png(paste(data.type,'_MDS_normalized.png', sep=""), width=960, height=768)
    plotMDS(d, col=as.numeric(d$samples$group), main=paste("MDS Plot normalized on", title.text, sep=" "))
    legend("bottomleft", as.character(unique(d$samples$group)), col=1:nlevels(d$samples$group), pch=20)
    dev.off()
} else {
    png(paste(data.type,'_PCA_normalized.png', sep=""), width=1280, height=1024)
    plot(log(d$pseudo.counts))
    title("three samples needed for PCA. x-y plot generated")
    dev.off()
}
                    
#calculate prior.df for new version of edgeR, which no longer accepts prior.n
#residual.df = #samples - #groups
#prior.df = residual.df * prior.n
pdf = (nrow(targets) - nlevels(d$samples$group)) * as.numeric(smooth.val)

if (stat.type == "qCML" || data.type == "essentialgenes"){
    if (data.type == "essentialgenes"){
        disp.trend = "loess"
    } else {
        disp.trend = "none"
    }
    disp.trend 
    if (disp.type == "tagwise"){
        if (nrow(targets) >2){
            d <- estimateTagwiseDisp(d, prior.df=pdf, trend=disp.trend)
            de.tagwise <-exactTest(d, dispersion="tagwise")
        } else {
            de.tagwise <-exactTest(d, dispersion="common")
        }
    }
    if (disp.type == "common") de.tagwise <-exactTest(d, dispersion="common")
} else {
    design <- model.matrix(~ group, data = targets)
    colnames(design) <- levels(d$samples$group) ## probably not necessary    
    d <- estimateDisp(d, design, prior.df=pdf, method="trend")
    if (disp.type == "common")  glmfit.tgw <- glmFit(d, design, dispersion = d$trended.dispersion)
    if (disp.type == "tagwise") glmfit.tgw <- glmFit(d, design, dispersion = d$tagwise.dispersion)
    de.tagwise <- glmLRT(glmfit.tgw)
}

toptags_tgw.out <-topTags(de.tagwise, n=nrow(de.tagwise)+1, adjust.method=adj.type)

## Generate volcano plot 
library(EnhancedVolcano)
png(paste(data.type,'_volcano.png', sep=""), width=1280, height=1024)
EnhancedVolcano(toptags_tgw.out$table, 
    lab=rownames(toptags_tgw.out$table), 
    x='logFC', 
    y='FDR',
    title="Volcano Plot",
    subtitle=data.type,
    ylab=bquote(~-Log[10]~ 'FDR'), 
    drawConnectors = T,
    legendPosition="bottom")
dev.off()

write.table(d$samples, paste(data.type, '_samples.tsv', sep=""), sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)
detags200 <- rownames(topTags(de.tagwise, n = 200)$table)
png(paste(data.type, '_plot.png', sep=""), width=1280, height=1024)
plotSmear(d, de.tags = detags200, main = paste('Fold change vs signal plot on', title.text, sep=" "))
abline(h = c(-2, 2), col = 'dodgerblue')
dev.off()
d$counts <- d$counts-1

samplenames <- read.delim("samplenames.txt", header=F)
samplenumber <- 1
while (samplenumber <= nrow(samplenames)) {
    colnames(d$counts)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
	colnames(d$pseudo.counts)[samplenumber] <-  as.matrix(samplenames)[samplenumber]
    samplenumber <- samplenumber + 1
}
if (data.type == "essentialgenes"){
    colnames(d$counts)[(nrow(samplenames)+1)] <-  "unique_flanking_sequences"
    colnames(d$pseudo.counts)[(nrow(samplenames)+1)] <-  "expected_reads"
}

ptt.file <- "genome.ptt"
if (data.type == "ta") ptt.file <- "ta.ptt"
allcounts.merged <- merge(d$counts, round(d$pseudo.counts), by=0)
row.names(allcounts.merged) <- allcounts.merged[, 1]
allcounts.merged <- allcounts.merged[, 2:(2 * nrow(targets)+1)] 
alldata.merged <- merge(allcounts.merged, toptags_tgw.out$table, by=0)
genefunctions <- read.delim(ptt.file, header=T) 
alloutput.merged <- merge (alldata.merged, genefunctions, by=1, all=TRUE)
if (data.type == "ta"){
    alloutput.merged <- alloutput.merged[, c(1,1,2:ncol(alloutput.merged))]
    alloutput.merged$Row.names.1 <- abs(as.numeric(alloutput.merged$Row.names.1))
    colnames(alloutput.merged)[2] <- "Position"
}
write.table(alloutput.merged, paste(data.type, '_alloutputmerged.tsv', sep=""), sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE) 
 
#density plot and possible cutoff calculation
library(zoo)

nr.x.minima <- 100
densityadjust <- 0.8 ## gene and ta
if (data.type == "essentialgenes") densityadjust <- 0.1
max.minima <- 6
if (data.type == "essentialgenes") max.minima <- 1
if (data.type == "ta") max.minima <- 4

while(nr.x.minima > max.minima) {
    d.density<-density(toptags_tgw.out$table$logFC, kernel = c("gaussian"), bw="nrd", n=2048, adjust=densityadjust)
    y <- d.density$y
    x <- d.density$x
    xz <- as.zoo(y)
    rxzmin <- rollapply(xz, 3, function(x) which.min(x)==2) #local minima
    rxzmax <- rollapply(xz, 3, function(x) which.max(x)==2) #local maxima (don't need)
    x.minima<- x[index(rxzmin)[coredata(rxzmin)]]
    y.minima<- y[index(rxzmin)[coredata(rxzmin)]]
    nr.x.minima=nrow(as.matrix(x.minima))
    densityadjust <- densityadjust + 0.1 
}
if (data.type == “essentialgenes”) write.table(round(x.minima,2),“logFC.minima.txt”,row.names = F, col.names = F)

y.range <- max(y)-min(y)
y.shift <- y.range/20
png(paste(data.type, '_densityplot.png', sep=""), width=1280, height=1024)
plot(d.density, main="Density plot of Log Fold change. Red dot depicts putative fold change cutoff(s)", xlab="Log2 fold change")
points(x.minima,y.minima, col="red", pch=19:21 )
text(x.minima,y.minima+y.shift, labels=round(x.minima,2))
dev.off()
q()
