Args <- commandArgs(TRUE)
counts <- read.delim(Args[1], header=F)
counts[(3)] <-abs(counts[(2)])
attach(counts)
sort.counts <-counts[order(V2.1) , ]
detach(counts)
y <- sort.counts[, (1)]
x <- sort.counts[, (3)]

## Remove extreme outliers prior to loess transformation 
coeff <- 3
Q <- quantile(y, probs=c(.25, .75), na.rm=F)
iqr <- IQR(y)
dn <- Q[1] - iqr*coeff
up <- Q[2] + iqr*coeff
outliers <- !(y<=up & y>=dn)
## Hack to make sure we always predict the full range of values
outliers[1] <- FALSE
outliers[length(outliers)] <- FALSE

y.filt <- y[!outliers]
x.filt <- x[!outliers]
y.loess <- loess(y ~ x, span=1, data.frame(x=x.filt, y=y.filt), control = loess.control(statistics = c("approximate"),trace.hat = c("approximate")))
y.predict <- predict(y.loess, data.frame(x=x))
y.ratio <- y.predict/median(y.predict)
y.new <- y/y.ratio
newcounts<- as.data.frame(y.new)
newcounts[(2)] <- sort.counts[(2)]
write.table(newcounts, Args[1], sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

## Cube root transform (to maintain zero and negative values) and plot
y.trans <- sign(y) * (abs(y) ^ (1/3))
y.predict.trans <- sign(y.predict) * (abs(y.predict) ^ (1/3))
y.new.trans <- sign(y.new) * (abs(y.new) ^ (1/3))
png(paste(Args[1], '_loess.png', sep=""), width=1280, height=1024)
plot(y.trans ~ x, pch=20, main=Args[1], xlab="Site", ylab="Counts ^ (1/3)", ylim=range(y.trans, y.predict.trans, y.new.trans))
points(x, y.new.trans, pch=4, col="red")
lines(x, y.predict.trans, col="blue")
dev.off()
