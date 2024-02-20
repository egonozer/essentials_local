Args <- commandArgs(TRUE)
genecounts <- read.delim(Args[1], header=T)
genelocations <- read.delim("genome.ptt", header=T)
row.names(genelocations)<-(genelocations[,(1)])
row.names(genecounts)<- (genecounts[,(2)])
fcandlocations <- merge (genecounts, genelocations, by=0)
attach(fcandlocations)
sort1.fcandlocations <- fcandlocations[order(Start) , ]
detach(fcandlocations)
x <-sort1.fcandlocations[,(5)]
y <-sort1.fcandlocations[,(2)]

## Remove extreme outliers prior to loess transformation 
coeff <- 3
Q <- quantile(y, probs=c(.25, .75), na.rm=F)
iqr <- IQR(y)
dn <- Q[1] - iqr*coeff
## Make sure 0 values are considered outliers since these are likely essential genes and shouldn't be included in loess calculation
## 0's also likely shouldn't be loess corrected, though the efffect of this is likely negligible
dn <- max(dn, 0)
up <- Q[2] + iqr*coeff
outliers <- !(y<=up & y>=dn)
## Hack to make sure the full range of values is predicted
outliers[1] <- FALSE
outliers[length(outliers)] <- FALSE

y.filt <- y[!outliers]
x.filt <- x[!outliers]
y.loess <- loess(y ~ x, span=1, data.frame(x=x.filt, y=y.filt))
y.predict <- predict(y.loess, data.frame(x=x))
y.ratio <- y.predict/median(y.predict)
y.new <- y/y.ratio
newcounts<- as.data.frame(y.new)
newcounts[(2)] <- sort1.fcandlocations[(1)]
colnames(newcounts) <- colnames(genecounts)
write.table(newcounts, Args[1], sep = "\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

## Cube root transform (to maintain zero values) and plot
y.trans <- sign(y) * (abs(y) ^ (1/3))
y.predict.trans <- sign(y.predict) * (abs(y.predict) ^ (1/3))
y.new.trans <- sign(y.new) * (abs(y.new) ^ (1/3))
png(paste(Args[1], '_loess.png', sep=""), width=1280, height=1024)
plot(y.trans ~ x, pch=20, main=Args[1], xlab="Gene", ylab="Counts ^ (1/3)", ylim=range(y.trans, y.predict.trans, y.new.trans))
points(x, y.new.trans, pch=4, col="red")
lines(x, y.predict.trans, col="blue")
dev.off()