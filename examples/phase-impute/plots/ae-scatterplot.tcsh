#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set baseind = ../results/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all.renamed-cols
set baseoutd = ../results.plots/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all

mkdir -p $baseoutd

# adapted from https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
# do over all individuals; tag counts must be 5+ on each haplotype not to have *** in the name
set inp = $baseind/final.ae.best.with-header.uniqs.txt
set tmp = $baseoutd/$inp:t.tmp
awk '$1 !~ /^\*/' $inp \
 >! $tmp
set inp = $tmp

# do individuals independently
R --no-save > /dev/stderr << __R__
  phaser = read.delim("$inp")
  ids <- unique(phaser[['sample.id']])
  for (id in ids) {
    print(id)
    pidv <- phaser[phaser[['sample.id']] == id,]
    phaserbinomp = apply(pidv[,c("refCount","altCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)[["p.value"]])
    phaserbinomq = p.adjust(phaserbinomp, method = "fdr")
    s <- cbind(pidv, phaserbinomq)
    s <- s[s[["phaserbinomq"]]<0.01,]
    write.table(s, file=paste("$baseoutd/", id, ".fdr0.01.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

    # plot haplotype A versus haplotype, highlight sites with significant imbalance (FDR < 5%)
    png(paste("$baseoutd/", id, ".scatterplot.png", sep=""))
    plot(pidv[["altCount"]], pidv[["refCount"]], log="xy", col=(phaserbinomq<0.01)+1, ylab="Haplotype B Count", xlab="Haplotype A Count", pch='.', main=id)
    abline(0,1,col="grey")
    a <- phaserbinomq[phaserbinomq<0.01]
    sz1 <- length(phaserbinomq) - length(a)
    sz2 <- length(a)
    legend("topleft",c(paste("No significant imbalance:", sz1), paste("Significant imbalance (FDR 1%):", sz2)),pch=c(15,15),col=c(1,2))
    dev.off()
  }
__R__

rm -f $tmp

exit 0
