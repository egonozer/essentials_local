Args <- commandArgs(TRUE)
counts <- read.delim(Args[1], header=F)
counts[(3)] <-abs(counts[(2)])
attach(counts)
sort.counts <-counts[order(V2.1) , ]
detach(counts)
y <- sort.counts[, (1)]
x <- sort.counts[, (3)]
#y.loess <- loess(y ~ x, span=1, data.frame(x=x, y=y))
y.loess <- loess(y ~ x, span=1, data.frame(x=x, y=y), control = loess.control(statistics = c("approximate"),trace.hat = c("approximate")))
y.predict <- predict(y.loess, data.frame(x=x))
y.ratio <- y.predict/median(y.predict)
y.new <- y/y.ratio
newcounts<- as.data.frame(y.new)
newcounts[(2)] <- sort.counts[(2)]
write.table(newcounts, Args[1], sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
