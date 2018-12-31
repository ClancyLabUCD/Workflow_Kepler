#!/usr/bin/python
import sys
#MSWF Reporting file parser

def f(key):
    return True


inputFile=sys.argv[1] # input file path
stageName=inputFile.split("/") # conf file name from input file path
with open(inputFile) as myConfFile: # open file to Read
    for line in myConfFile:
	if "=" in line:
		confName, confValue = line.partition("=")[::2]         
		print(confName+' = '+confValue)

