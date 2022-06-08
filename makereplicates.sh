#!/bin/bash

path1=$1

desiredcopies=${2:-3} # number of desired copies. Defaults to three replicates = two copies
current_copy=1

while [ $current_copy -ne $desiredcopies ]
do
    current_copy=$(( $current_copy +1 ))
    targetdir=${path1::-2}$current_copy
    echo "copying from ${path1} to ${targetdir}"
    cp -rf ${path1} ${targetdir}

    
    echo "copy ${current_copy} done"
done
echo "done making replicates"