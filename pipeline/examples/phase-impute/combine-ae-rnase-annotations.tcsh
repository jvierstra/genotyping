#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set dstamp = 2019-04-04
set baseind = ../results/$dstamp/rnase-ae/output/expression

foreach vcfd (`find $baseind/ -maxdepth 1 -mindepth 1 -type d`)
  set fs = ()
  foreach sample (`find $vcfd/ -maxdepth 1 -mindepth 1 -type d`)
    set chroms = (`find $sample/ -maxdepth 1 -mindepth 1 -type f -name "*.ae.txt"`)
    if ( $#chroms > 0 ) then
      set outd = $chroms[1]:h
      set header = $outd/.header
      awk 'NR == 1' $chroms[1] \
       >! $header

      cat $chroms \
        | awk '$2 != "start"' \
        | sort-bed - \
        | cat $header - \
       >! $outd/final.ae.with-header.bed

      chmod 755 $outd/final.ae.with-header.bed

      rm -f $header

      set fs = ($fs $outd/final.ae.with-header.bed)
    endif
  end

  if ( $#fs > 0 ) then
    set outd = $fs[1]:h
    set header = $outd/.header
    awk 'NR == 1' $fs[1] \
     >! $header

    cat $fs \
      | awk '$2 != "start"' \
      | sort-bed - \
      | cat $header - \
     >! $vcfd/final.ae.all-samples.with-header.bed

    awk '$8 !~ /inf/' $vcfd/final.ae.all-samples.with-header.bed \
      | awk '$8 != 0.0' \
     >! $vcfd/final.ae.best.with-header.bed

    awk 'NR > 1' $vcfd/final.ae.best.with-header.bed \
      | sort -gk8,8 --buffer-size=16G \
      | (head -500; tail -500) \
      | cat $header - \
     >! $vcfd/final.ae.tophits.with-header.bed

    chmod 755 $vcfd/final.ae.all-samples.with-header.bed
    chmod 755 $vcfd/final.ae.best.with-header.bed

    rm -f $header
  endif
 
end

exit 0
