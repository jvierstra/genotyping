#!/bin/bash -x

cd "$(dirname "$0")"
tag=genotyping-100-donors-t-cells-selected-hq-hg19
datadir=../data/input.no36bpfilter

mkdir -p $datadir

source "/home/solexa/stampipes/scripts/lims/api_functions.sh"

data=$(lims_get_all "aggregation/file_detail/?tag_slug=$tag&page_size=1000")

jq -r '[ "AG" + (.id | tostring), .files["all-alignments-bam"], .files["hotspot-calls"]] | join(" ")' <<< "$data" \
  > files.txt

rm -rf $datadir/*bam $datadir/*bai

while read ag bam spots ; do
	echo "$ag" "$bam" "$spots"
	ln -s "$bam" "$datadir/$ag.bam"
	ln -s "$bam.bai" "$datadir/$ag.bam.bai"
	#ln -s "$spots" "$datadir/$ag.spots.starch"
done < files.txt

exit 0
