#!/usr/bin/env bash

## Usage: bash csv2vcf.sh input.csv output.vcf
## Converts Varvis candidate variants (.csv) to vcf
## As of today (2020-11-27) Varvis export uses hg19!!!

upload_id=$(basename ${1%.tsv})
upload_id_fix=$(echo $upload_id | tr "_" "-")
##person_id=$(echo $upload_id | tr "_" "\t" | cut -f1)
##echo $person_id

##tempdir="$(dirname $2)/temp"
##[ -d $tempdir ] || mkdir -p $tempdir

tempdir=$(mktemp -d)
echo "TEMPDIR: $tempdir"

## csv to tsv
## cat $1 | dos2unix | sed 's/,/\t/g' | sed 's/^\t/NA\t/; :a;s/\t\(\t\|$\)/\tNA\1/;ta' > $tempdir/${upload_id}.tsv

## map
tail -n+2 "$1" \
    | awk 'BEGIN{OFS="\t"}{id=sprintf("%s_%s_%s_%s_%s", $1, $2, $3, $4, $5); print $1, id, 0, $2}' \
    | sed 's/^XY/25/' \
    | sed 's/^MT/26/' \
    | sed 's/^X/23/' \
    | sed 's/Y/24/' \
    > $tempdir/${upload_id}.map

### ped
tail -n+2 "$1" | awk -v upload_id_fix=$upload_id_fix 'BEGIN{OFS=" "; aggr=""}{aggr=aggr " " $4 " " $5}END{print upload_id_fix, upload_id_fix, 0, 0, 0, -9 aggr}' > $tempdir/${upload_id}.ped

### convert map + ped to vcf
plink --file $tempdir/$upload_id --recode vcf --out $tempdir/$upload_id

cat $tempdir/${upload_id}.vcf \
    | sed 's/^23/X/' \
    | sed 's/^24/Y/' \
    | sed 's/^25/XY/' \
    | sed 's/^26/MT/' \
    | sed 's/##contig=<ID=23,/##contig=<ID=X,/' \
    | sed 's/##contig=<ID=24,/##contig=<ID=Y,/' \
    | sed 's/##contig=<ID=25,/##contig=<ID=XY,/' \
    | sed 's/##contig=<ID=26,/##contig=<ID=MT,/' \
    > $tempdir/${upload_id}.temp.vcf

mv $tempdir/${upload_id}.temp.vcf "$2"

rm -r $tempdir
