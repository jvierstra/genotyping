#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set baseind = ../results/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all.renamed-cols
set baseoutd = ../results.plots/$dstamp/ai-calc-rnaseq-ae.filteredR2/output/expression/filtered.all

mkdir -p $baseoutd

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
    phaserdistr = apply(pidv[,c("refCount","altCount")], 1, function(x) x[2]/(x[1]+x[2]))

    png(paste("$baseoutd/", id, ".distr.png", sep=""))
    hist(phaserdistr, col="peachpuff", border="black", prob=T, main=id, xlab="ref/(ref+alt)", axes=F, xlim=c(0,1), ylab="")
    lines(density(phaserdistr), lwd=2, col="cornflowerblue")
    axis(1)
    dev.off()
  }
__R__

rm -f $tmp

exit 0
