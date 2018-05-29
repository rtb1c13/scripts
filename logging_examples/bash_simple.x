#!/bin/bash

# Script for various simple argument parsing examples
# For this example, first line should be the path to a config file
# that defines variables $foofromfile and $barfromfile

# Positional arguments
foo=$1
bar=$2
echo "Argument 1 is $foo" > logfile_simpleargs.txt
echo "Argument 2 is $bar" >> logfile_simpleargs.txt


# Source argument 1 - it can be a config file!
. $foo
echo "First variable in file $foo is $foofromfile" >> logfile_simpleargs.txt
echo "Second variable in file $foo is $barfromfile" >> logfile_simpleargs.txt

# Or read lines and perform actions
while read line
do
 echo "Now reading line $line from file $foo" >> logfile_simpleargs.txt
done < $foo
