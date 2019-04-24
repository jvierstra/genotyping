#!/bin/tcsh -efx
# author : sjn
# date : Apr.2019

set csvs = (input/production_qc_dnase_20190424-120705.csv input/production_qc_rna_20190424-120657.csv)

foreach csv ($csvs)
  cat $csv \
    | dos2unix \
    | sed 's; ;_;g' \
    | awk -F"," \
        '{ \
           for(i=1;i<=NF;++i) { \
             if ($i ~ /\"/) { \
               if(cnt%2==1) { printf ","; } \
               cnt+=1; \
             } else { \
               if(i>1) { printf ","; } \
               printf "%s", $i; \
             } \
           } \
           printf "\n" \
         }' \
    >! $csv.tmp

 set allempty_spaces = \
                `awk -F"," \
                  '{ \
                    for(i=1;i<=NF;++i) { \
                      if ( $i == "" ) {empty[i]+=1} \
                    } \
                  } END { \
                    for(e in empty) { \
                      if(empty[e]==NR-1) { \
                        print e; \
                      } \
                    } \
                  }' $csv.tmp \
                 | sort -n \
                 | uniq \
                 | tr '\n' ',' \
                 | sed -r 's/(.*),/\1/'`

 set anyempty_spaces = \
              `awk -F"," \
                  '{ \
                    for(i=1;i<=NF;++i) { \
                      if ( $i == "" ) { print i; } \
                    } \
                  }' $csv.tmp \
                 | sort -n \
                 | uniq \
                 | tr '\n' ',' \
                 | sed -r 's/(.*),/\1/'`

  # if you want to remove columns with quotes
  #  and those that are empty across all rows
#  cut -f$allempty_spaces --complement -d"," $csv.tmp \
#   >! $csv:r.parsed.csv

  # if you want to remove columns with quotes
  #  and those that are empty in any row
  cut -f$anyempty_spaces --complement -d"," $csv.tmp \
   >! $csv:r.parsed2.csv

  rm $csv.tmp
end

exit 0
