#!/bin/bash




Help()
{
    echo
    echo "Transformato helper scripts suite: tf_makereplicates"
    echo "Script to make copy an arbitrary number of replicate folders"
    echo
    echo "-------------USAGE-------------------------------"
    echo "USAGE: tf_analyse [PATH_TO_REPLICATE_FOLDER] [DESIRED_NUMBER_OF_REPLICATES]"
    echo
    echo "-h    Prints this help"
    echo 
    echo
    echo "By default, two copies are created. [DESIRED_NUMBER_OF_REPLICATES] specifies the total amount of folders you want"
    echo "So if you specify you want 6 replicates, it wll make 5 copies: -2, -3, -4, -5, -6"
    echo
    echo "------------FILE HANDLING---------------------------"
    echo "Source folder must be specified as relative path, ending with -1/"
    echo "Examples of correct commands:"
    echo "tf_makereplicates ./replicate-1/ #Makes 2 replicates (replicate-2 and replicate-3) in the local directory"
    echo "tf_makereplicates ~/tf_production/struc22/23to28/5ns_k3ex3-scaling-1/ 10 # Makes 10 replicates (-2 to -10) in the same location as the source"
}

while getopts ":h" option; do
    case $option in
        h) # Help
            Help
            exit;;
        
    esac
done

path1=$1

desiredcopies=${2:-3} # number of desired copies. Defaults to three replicates = two copies
current_copy=1

if [[ ${path1: -1} != "/" ]]
then
    echo "No trailing slash - adding"
    
    path1="${path1}/"

fi

while [ $current_copy -ne $desiredcopies ]
do
    current_copy=$(( $current_copy +1 ))
    targetdir=${path1::-2}$current_copy
    echo "copying from ${path1} to ${targetdir}/"
    cp -rf ${path1} ${targetdir}

    
    echo "copy ${current_copy} done"
done
echo "done making replicates"