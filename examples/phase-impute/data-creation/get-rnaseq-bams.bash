#!/bin/bash -e
# author : sjn
# date : Apr.2019

if [ $# -ne 2 ]
then
  printf "Wrong number of args\n"
  exit -1
fi

tagid=$1
output=$2

source ~sjn/bin/getlims.bash

vs="jq -r '. | { \
      allbam: .files.all_alignments_bam \
    }'"

# jq does not like '-' in keys when filtering ie; .stats.nuclear-align
lims_get_all "aggregation/file_detail/?tag_slug=$tagid" \
  | awk -F":" \
    'BEGIN {OFS=FS} ; { \
      gsub(/-/, "_", $1); print; \
    }' \
  | eval $vs \
  | awk '{ if ($1 == "\"allbam\":") { print $2; } }' \
  | sed 's;";;g' \
 >$output

exit 0
