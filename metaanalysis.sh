#!/bin/bash

metafile="metaanalysis.txt"
metacsv="metaanalysis.csv"


Help()
{
    echo "Transformato helper scripts suite: tf_analysis"
    echo "Script to retrieve run data from analysis_*.out files"
    echo
    echo "-------------USAGE-------------------------------"
    echo "USAGE: tf_analyse [-i CSVNAME|-h]"
    echo
    echo "-h    Prints this help"
    echo "-o    Redirects the output csv form metaanalysis.csv to your name"
    echo
    echo "Run parameters are retrieved via regex from filenames."
    echo "Relative free energies are retrieved from files themselves"
    echo
    echo "------------FILE HANDLING---------------------------"
    echo "Files must be in current working directory"
    echo "To be recognized, filenames must have the following form:"
    echo "analysis_[TIMESCALE]-[KVALUE][RESTRAINTTYPE]-[REPLICATE].out"
    echo "Example:"
    echo "analysis_10ns-k100ex3-2.out"
}

while getopts ":h:o:" option; do
    case $option in
        h) # Help
            Help
            exit;;
        o)
            metacsv=${OPTARG};;
        \?)
            echo "Unrecognized parameter. Showing help instead."
            Help
            exit 1;;
    esac
done

echo Fetching Metaanalysis data: All of analysis_**.out and writing to "${metafile}"

echo "Meta - Analysis" > "${metafile}"

for i in analysis_**.out;
do echo ${i}>>"${metafile}";
cat ${i} | grep uncertanty >> "${metafile}";
done;

echo  Fetching Metaanalysis data: All of analsis_**.out and creating csv at "${metacsv}"


echo "name,timescale,kvalue,restraint,ktvalue,kcalvalue"> ${metacsv}

for i in analysis_**.out;
do ktvalue=$(perl -ne 'm/(?:(?!.*\+-.*)^Free energy difference: (.*) \[kT\]$)/ && print "$1\n"' ${i}); # fucking perl becaus grep and sed are piles of shit
timescale=$(echo ${i} | perl -ne '/(?:analysis_([\d\.]*ns))/g && print "$1\n"')
kvalue=$(echo ${i} | perl -ne '/(?:.*_.*-(k.*)(?:\Qex\E|\Qsimple\E))/g && print "$1\n"')
restraint=$(echo ${i} | perl -ne '/(?:.*_.*-.*(\Qex\E\d*|\Qsimple\E|\Qnorestraint\E))/g && print "$1\n"')


addline="${i},${timescale},${kvalue},${restraint},${ktvalue}"
echo ${addline}
echo ${addline} >> ${metacsv};
done;


echo "done";

