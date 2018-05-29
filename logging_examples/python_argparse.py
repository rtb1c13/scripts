#/usr/bin/env python

# Argparse & from file example

import sys, argparse


### Argparser
def parse():
   # fromfile_prefix_chars allows params to be read from file
   parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
   parser.add_argument("-s","--sequence",help="Sequence for analysis", nargs='+',required=True)
   parser.add_argument("-n","--number",help="Number of repeats for sequence", required=True, type=int)
   parser.add_argument("-l","--logfile",help="Logfile. Default: argparse_logs.txt", type=str, default='argparse_logs.txt')
   args = parser.parse_args()
   return args

### Now the main script
arguments = parse()

# argparse won't split multiple arguments on a line - you might need something like the below.
# Alternatively create a special version of parser.convert_arg_line_to_args in the parse function.
if len(arguments.sequence) == 1:
    arguments.sequence = arguments.sequence[0].split()

output = arguments.sequence * arguments.number
with open(arguments.logfile, 'a') as f:
    f.write("Command-line arguments: "+' '.join(i for i in sys.argv)+'\n')
    f.write("Output : %s\n" % output)
