#!/bin/bash 
if [ -z $3 ]
then
    echo "USAGE: bash $0 1_dataDir 2_suffix 3_numberOfMinGenes\n i.e., bash $0 data xls 500 \n"
    exit
fi

if [ ! -d $1 ]
then 
    echo "ERROR: $1 doesn't exist."
    exit
fi

fs=`ls $1/*$2`
echo $fs

for  f in $fs 
do
    bn=`basename $f`
    ofPrefix=${bn/%${2}/}
    tdir=${1}/$ofPrefix
    if [ -d $tdir ]
    then
	continue
    fi
    echo "Rscript ~/codes/git/SingleCell/basic_analysis.R $f $ofPrefix $3"
    Rscript ~/codes/git/SingleCell/basic_analysis.R $f $ofPrefix $3
done


exit
