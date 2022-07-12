#!/bin/bash

READ1=$1
READ2=$2
REF=$3
OUTPUT=$4
THREADS=$5
SAMPLE=$6

echo "In .sh file";
ID=$(gunzip -c "$READ1" | head -n 1 | cut -c 2- | awk -v FS=':' '{print "ID:"$1"."$2}');
echo $ID;
echo post_id;
PU=$(gunzip -c "$READ1" | head -n 1 | cut -c 2- | awk -v FS=':' '{print "PU:"$3"."$4}');
echo $PU;
echo post_pu;

echo "$SAMPLE"
echo '@RG\tSM:'$SAMPLE'\tLB:'$SAMPLE'\t'$ID'\t'$PU'.'$SAMPLE'\tPL:ILLUMINA'

bwa mem $REF $READ1 $READ2 \
    -t $THREADS \
    -R '@RG\tSM:'$SAMPLE'\tLB:'$SAMPLE'\t'$ID'\t'$PU'.'$SAMPLE'\tPL:ILLUMINA' \
| samtools view -Sb - > $OUTPUT