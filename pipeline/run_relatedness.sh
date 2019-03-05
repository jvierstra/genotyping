#!/bin/bash -x

usage () {
  echo -e "Usage: $0 <input.vcf.gz> <relative-thold> <outdir>" >&2
  exit 2
}

if [[ $# != 3 ]]; then
   usage
fi

input=$1
relative_thold=$2
outdir=$3

mkdir -p $outdir

source /net/module/Modules/default/bash
module load vcftools

output=$outdir/$input:t:r.relatedness
vcftools --gzvcf $input --relatedness --stdout \
  >! $outdir/$input:t:r.relatedness

# clean up and merge related individuals; produce output based upon transitive closure and an
#  alternative based upon greedy selections

# Transitive closure.  If A-B and B-C, then A-C even if that relationship does not meet threshold
awk '$NF != "-nan"' $output \
  | awk '$1 != $2' \
  | awk -v t=$relative_thold \
      'BEGIN {OFS="\t"; a[-1]; r[-1]; all[-1]; delete a[-1]; delete r[-1]; delete all[-1]} ; { \
        if ( NR > 1 ) { \
          all[$1]; all[$2]; \
          if ( $NF > t ) { \
            if ( $1 in r ) { \
              a[r[$1]] = a[r[$1]]";"$2; \
              r[$2] = r[$1]; \
              delete a[$1]; \
              if ( $2 in a ) { \
                a[r[$1]] = a[r[$1]]";"a[$2]; \
              } \
              delete a[$2]; \
            } else { \
              r[$2] = $1; \
              if ( $1 in a ) { \
                a[$1] = a[$1]";"$2; \
              } else { \
                a[$1] = $2; \
              } \
              if ( $2 in a ) { \
                a[$1] = a[$1]";"a[$2]; \
              } \
              delete a[$2]; \
            } \
          } \
        } \
      } END { \
        i=1; \
        for ( j in all ) { \
          if ( j in a ) { \
            print j, i; \
            n=split(a[j], b, ";"); \
            uniqs[j]; \
            for(k=1;k<=n;++k) { \
              if ( !(b[k] in uniqs) ) { \
                print b[k], i; \
			      uniqs[b[k]]; \
              } \
            } \
            ++i; \
          } else if ( !(j in r) ) { \
            print j, i; \
            uniqs[j]; \
            ++i; \
          } \
        } \
      }' \
 > $output.tc

# first come (by significance), first serve
awk 'NR > 1' $output \
  | awk '$NF != "-nan"' \
  | awk '$1 != $2' \
  | sort -grk3,3 \
  | awk -v t=$relative_thold \
      'BEGIN {OFS="\t"; a[-1]; take[-1]; delete a[-1]; delete take[-1]} ; { \
        if ( $NF > t ) { \
          if ( !($1 in take) && !($2 in take) ) { \
            if ( $1 in a ) { \
              if ( $2 in a ) { \
                a[$1] = a[$1]";"a[$2]; \
                delete a[$2]; \
              } else { \
                a[$1] = a[$1]";"$2; \
              } \
              take[$2]; \
            } else if ( $2 in a ) { \
              a[$2] = a[$2]";"$1; \
              take[$1]; \
            } else { \
              a[$1] = $1";"$2; \
              take[$2]; \
            } \
          } \
        } else { \
          if ( !($1 in take) ) { \
            if ( !($1 in a) ) { \
              a[$1] = $1; \
              take[$1]; \
            } \
          } \
          if ( !($2 in take) ) { \
            if ( !($2 in a) ) { \
              a[$2] = $2; \
              take[$2]; \
            } \
          } \
        } \
      } END { \
        i=1; \
        for ( j in a ) { \
          n=split(a[j], b, ";"); \
          print j, i; \
          for(k=1;k<=n;++k) { \
            if ( j != b[k] ) { \
              print b[k], i; \
            } \
          } \
          ++i; \
        } \
      }' \
 > $output.greedy

exit 0
