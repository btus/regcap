#!/bin/sh

#Makes hours, days and months directories and moves all files with those extensions
#to the appropriate directories. 
#Enter 'sh <file_path>/TimeSeriesMove.sh <path_to_target_directory>' in the terminal. Done. 
#For Example: sh /Users/brennanless/GoogleDrive/regcap/R_utilities/TimeSeriesMove.sh /Users/brennanless/GoogleDrive/regcap/out/humidity_controls_wDH/All_Combinations/

cd $1
mkdir hours
mkdir days
mkdir months
mv *.hours hours
mv *.days days
mv *.months months