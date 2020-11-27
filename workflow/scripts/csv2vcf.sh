#!/usr/bin/env bash

## Usage: bash csv2vcf.sh input.csv output.vcf
## Converts Varvis candidate variants (.csv) to vcf
## As of today (2020-11-27) Varvis export uses hg19!!!

bname=$(basename ${1%.csv})
person_id=$(echo $bname | tr "_" "\t" | cut -f1)

##tempdir="$(dirname $2)/temp"
##[ -d $tempdir ] || mkdir -p $tempdir

tempdir=$(mktemp -d)
echo "TEMPDIR: $tempdir"

## csv to tsv
cat $1 | dos2unix | sed 's/,/\t/g' | sed 's/^\t/NA\t/; :a;s/\t\(\t\|$\)/\tNA\1/;ta' > $tempdir/${bname}.tsv 

## map
tail -n+2 $tempdir/${bname}.tsv \
    | awk 'BEGIN{OFS="\t"}{id=sprintf("%s_%s_%s_%s_%s", $3, $4, $5, $6, $7); print $3, id, 0, $4}' \
    | sed 's/^XY/25/' \
    | sed 's/^MT/26/' \
    | sed 's/^X/23/' \
    | sed 's/Y/24/' \
    > $tempdir/${bname}.map

## ped
tail -n+2 $tempdir/${bname}.tsv | awk -v person_id=$person_id 'BEGIN{OFS=" "; aggr=""}{aggr=aggr " " $6 " " $7}END{print person_id, person_id, 0, 0, 0, -9 aggr}' > $tempdir/${bname}.ped

# convert map + ped to vcf
plink --file $tempdir/$bname --recode vcf --out $tempdir/$bname

cat $tempdir/${bname}.vcf \
    | sed 's/^23/X/' \
    | sed 's/^24/Y/' \
    | sed 's/^25/XY/' \
    | sed 's/^26/MT/' \
    | sed 's/##contig=<ID=23,/##contig=<ID=X,/' \
    | sed 's/##contig=<ID=24,/##contig=<ID=Y,/' \
    | sed 's/##contig=<ID=25,/##contig=<ID=XY,/' \
    | sed 's/##contig=<ID=26,/##contig=<ID=MT,/' \
    > $tempdir/${bname}.temp.vcf

mv $tempdir/${bname}.temp.vcf $2

rm -r $tempdir
