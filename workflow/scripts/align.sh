#!/bin/bash

READ1=$1
READ2=$2
REF=$3
OUTPUT=$4
THREADS=$5
SAMPLE=$(echo $6 | cut -d _ -f 1)  # Take sample from run name
LIBRARY=$(echo $6 | cut -d _ -f 2)  # Take library from run name
FASTQ_FIRST_LINE=$7

FIRST_LINE=$(cat "$FASTQ_FIRST_LINE" | cut -c 2-);
ID=$(echo $FIRST_LINE | awk -v FS=':' '{print "ID:"$1"."$2}');
PU=$(echo $FIRST_LINE | awk -v FS=':' '{print "ID:"$3"."$4}');
echo $ID;
echo post_id;
echo $PU;
echo post_pu;
echo "$SAMPLE";
echo '@RG\tSM:'$SAMPLE'\tLB:'$SAMPLE'\t'$ID'\t'$PU'.'$SAMPLE'\tPL:ILLUMINA';

bwa mem $REF $READ1 $READ2 \
    -t $THREADS \
    -R '@RG\tSM:'$SAMPLE'\tLB:'$LIBRARY'\t'$ID'\t'$PU'.'$SAMPLE'\tPL:ILLUMINA' \
| samtools view -Sb - > $OUTPUT