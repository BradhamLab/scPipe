#! /bin/bash

# undo fastq segregation performed after fastp qc

#1 should be segregated directory

#2 should be output directory

# example ./undo_segregation.sh <segegrated_qc> <qc>

copy_fastq () {
    local input_dir=$1
    local output_dir=$2
    local n_samples=`ls $input_dir | wc -l`
    if [ $n_samples -gt 0 ];
    then 
        for sample in `find $input_dir/* -type d`
        do
            local id=`basename $sample`
            new_loc=$output_dir"/$id/"
            `cp -r $sample/*fastq.gz $new_loc`
        done
    fi
}

copy_fastq "$1/good" $2

copy_fastq "$1/ugly" $2

copy_fastq "$1/bad" $2
