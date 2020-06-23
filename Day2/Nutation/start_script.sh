#!/bin/bash

#Version 1
#This script executes the make_inputfile script and uses loops to scan the parameters given in the defined data input files.

if [ $# -eq 0 ]
  then
    echo "No arguments supplied. You have to specify a rf power."
fi

rf_power=$1


chmod +x make_inputfile.sh

for pulse_length in $(< pulse_length.in)
do
	./make_inputfile.sh nutation $rf_power $pulse_length
done

cat $(find ./ -name "*.xy" | sort -V) > nutationcurve.xy

gawk -i inplace '{print NR " " $2}' nutationcurve.xy