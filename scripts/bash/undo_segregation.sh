#! /bin/bash

# undo fastq segregation performed after fastp qc

#1 should be segregated directory

#2 should be output directory

# example ./undo_segregation.sh <segegrated_qc> <qc>


for sample in `find $1good/* -type d`
do
    id=`basename $sample`
    new_loc=$2"/$id/"
    `cp -r $sample/*.fastq.gz $new_loc`
done

for sample in `find $1bad/* -type d`
do
    id=`basename $sample`
    new_loc=$2"/$id/"
    `cp -r $sample/*.fastq.gz $new_loc`
done

for sample in `find $1ugly/* -type d`
do
    id=`basename $sample`
    new_loc=$2"/$id/"
    `cp -r $sample/*.fastq.gz $new_loc`
done