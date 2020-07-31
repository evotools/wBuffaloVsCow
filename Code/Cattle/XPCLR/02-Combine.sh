#!/bin/bash

infld=OUTPUTS
chrom=$1

mkdir XPCLR_RESULTS
for i in `ls ${infld}`; do 
    mkdir XPCLR_RESULTS/${i}; 
    for x in $( seq 1 $chrom ); do 
        if [ $x -eq 1 ]; then 
            cat ${infld}/$i/${i}.${x}.xpclr
        else 
            awk 'NR>1{print }' ${infld}/$i/${i}.${x}.xpclr
        fi
    done > XPCLR_RESULTS/${i}/${i}.xpclr 
    echo "Done ${i}"
done