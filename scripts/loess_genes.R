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
y.loess <- loess(y ~ x, span=1, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
y.ratio <- y.predict/median(y.predict)
y.new <- y/y.ratio
newcounts<- as.data.frame(y.new)
newcounts[(2)] <- sort1.fcandlocations[(1)]
colnames(newcounts) <- colnames(genecounts)
write.table(newcounts, Args[1], sep = "\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
