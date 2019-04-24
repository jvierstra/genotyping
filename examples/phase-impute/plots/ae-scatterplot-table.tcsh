#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set indir = results.plots/$dstamp/rnaseq-ae.filteredR2/output/expression/filtered.all
set blogdir = genotypes/genotyping/altius100donors/$indir
set b = (`find ../$indir/ -maxdepth 1 -mindepth 1 -type f -name "*.png" | xargs -I % sh -c 'basename % | cut -f1 -d"." | cut -f2 -d"-"'`)
set b = (`echo $b | tr ' ' '\n' | sort -g`)
set pngs = Indiv-{`echo $b | tr ' ' ','`}.scatterplot.png
set pngs = `echo $pngs | tr ' ' '\n' | awk '{ if (NR%6==0) { print s";"$1 } else if (NR%6==1) { s=$1; } else { s=s";"$1 } } END { print s }'`

set outdir = ../$indir

set html = $outdir/simple.html
printf '<html>\n<table border="1">\n' >! $html
foreach png_group ($pngs)
  printf '  <tr>\n' >> $html
  foreach png (`echo $png_group | tr ';' ' '`)
    set genes = $png:t:r:r.fdr0.01.txt
    printf '    <td><embed src="'$png'" /><a href="'$genes'">GeneList</a></td>\n' >> $html
  end
  printf '  </tr>\n' >> $html
end
printf '</table>\n</html>\n' >> $html

exit 0
