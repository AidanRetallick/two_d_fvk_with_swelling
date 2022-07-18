#!/bin/sh
# Quick script to batch the data aquisition
# The Directory
driver_name="circular_cap_fvk"
arguments="--p 10 --dp 10 --p_max 11"
areaprefix="0.0"
# Loop and run the code 
for i in `seq 1 6`;
do
    echo "##########################################"
    echo $circular_cap_fvk "--element_area" $areaprefix$i $arguments 
    ./$driver_name "--element_area" $areaprefix$i $arguments | tee RESLT/output$areaprefix$i.txt
done

# Now grab the data
grep -h "w in the middle: " RESLT/output0.* | grep -Eo "[0-9]\.[0-9]+" > RESLT/centres.dat
ls RESLT/output0.* | grep -Eo "[0-9]\.[0-9]+" > RESLT/areas.dat
# Plot centres versus areas in your favourite plotting package
