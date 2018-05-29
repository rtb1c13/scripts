#!/bin/bash

# Script for a getopts argument parsing example
# For this example, -a is a flag, -b defines a variable, and -h displays help

usage() { echo "Usage: $0 [-a] [-b value] -h for help"; exit 1; }

# First colon = silent error reporting.
if ( ! getopts ":ab:h" opt ); then
 usage
fi

# Write current command line to output file
echo $0 $@ >> getopts_logfile.txt

# Use 'case' switch to define behavior
while getopts ":ab:h" opt; do
 case $opt in
  a)
   echo "-a flag was triggered!" >> getopts_logfile.txt
   ;;
  b)
   echo "-b variable is $OPTARG" >> getopts_logfile.txt
   ;;
  h)
   echo "-h prints help"
   usage
   ;;
  \?)
   echo "Unknown option -$OPTARG" >> getopts_logfile.txt
   exit 1
   ;;
  :)
   echo "Option -$OPTARG expects an argument" >> getopts_logfile.txt
   exit 1
   ;;
 esac
done
