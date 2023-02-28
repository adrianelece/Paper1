#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Imports
#########

#from ast import pattern
import sys, os, glob
import numpy as np
from collections import defaultdict
import pandas as pd

# Defaults
##########

out_file  = "Data.tsv"
pattern_file = "network*.tsv"

__doc__ = """
SYNOPSIS

Python script to extract differentially methylated regions from individual samples

   DatabyDir.py -r <regions file> -s <directory to the samples> -p <pattern to found in the files> -o <output> [-h] [-v] [-a]

DESCRIPTION

 Parameters:
    -r,--region         txt file with all the regions in common. [Required]
    -s,--sample         Path to where the files are stored [Required]
    -o,--output         Output file name [Optional]
    -p, --pat           Pattern in the name of the files to found [Required]
    -v,--version        Shows version
    -h,--help           Shows help options
    -a,--author         Shows author information
 
"""
__author__ = "Adrian Lopez-Catalina"

# Functions
###########

def paramerror(param, message):
    """Gives an error message if argument is not valid"""
    print( "\tERROR: Parameter '", param, "' ", message, "\n")
    help()
    sys.exit(1)

def messagerror(message):
    """Gives an error message if something is wrong"""
    print( "\tERROR:", message, "\n")
    help()
    sys.exit(1)

def help():
    print( globals()['__doc__'])
    sys.exit(1)

def version():
    print( globals()['__version__'])
    sys.exit(1)

def author():
    print( globals()['__author__'])
    sys.exit(1)
    
# Parameters catching
#####################

if len(sys.argv) <= 1:
	messagerror("Parameters necessary. Please, use -h.")
else:
    script_name = sys.argv.pop(0)
    while len(sys.argv) > 0:
        param = sys.argv.pop(0)
        if param == '-o' or param == '--out':
            out_file = sys.argv.pop(0)
        elif param == '-r' or param == '--region':
            region_file = sys.argv.pop(0)
        elif param == '-p' or param == '--pat':
            pattern_file = sys.argv.pop(0)
        elif param == '-s' or param == '--sample':
            sample_file = sys.argv.pop(0)
        elif param == '-h' or param == '--help':
            help()
            sys.exit(1)
        elif param == '-a' or param == '--author':
            author()
            sys.exit(1)
        else:
            paramerror(param, "Invalid input")

if len(region_file) == 0:
    print("\tPlease, select your regions file.")
    help()
    sys.exit(1)
else:
    if len(sample_file) == 0:
        print("\tPlease, select your sample file")
        help()
        sys.exit(1)



path = sample_file
all_files = glob.glob(os.path.join(path , pattern_file))

test1=pd.read_csv(region_file, sep="\t",  header=None)
test1.columns=["region"]
allsam=pd.DataFrame(columns=["region","llr"])
dirc=''.join(path)
vcnames = []
for filename in all_files:
    colnametmp = filename   
    newcol = colnametmp.replace(dirc,"")
    vcnames.append(newcol)
    test2=pd.read_csv(filename, sep="\t",  header=None)
    test2.columns=["region","llr"]
    test = test1.merge(test2, how='left', on='region',indicator='merge')
    test.loc[test['merge'] == 'left_only', 'llr' ] = '.'
    out=test[["region","llr"]]
    allsam = allsam.merge(out,how='right', on='region')

cnamesempty =['na','region']
vcnamesfinal = [s.replace('.tsv',"") for s in vcnames]
cnamesempty=cnamesempty + vcnamesfinal
allsam.columns=cnamesempty
allsam.to_csv(out_file,index=False, sep="\t")
