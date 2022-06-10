#!/usr/bin/env python

#title	  :filter_CharGer.py
#author	  :Fernanda Martins Rodrigues (fernanda@wustl.edu)
#date	  :20220501

"""
	This script will filter CharGer output to contain only non-synonymous variants and variants with CharGer score ≥ 4. 
	It will also create a file containing only variants passing a given MAF threshold.

	Usage:
		python filter_CharGer.py -c [input CharGer file] -a [MAF threshold] -o [output file basename] -O [output directory]

	Arguments:
		-c, --charger:			charger output tsv file
		-a, --maf:	  			MAF threshold
		-O, --outputDirectory:  directory to write output files to

"""

import sys
import argparse
import getopt
import os
import re

def argument_parser():
	# create parser
	parser = argparse.ArgumentParser(description=__doc__)
	# add arguments
	parser.add_argument("-c", "--charger", required=True, help="charger output tsv file corresponding to input VCF file")
	parser.add_argument("-a", "--maf", required=True, help="MAF threshold")
	parser.add_argument("-o", "--outBasename", required=True, help="output file basename")
	parser.add_argument("-O", "--outputDirectory", required=True, help="directory to write output files to")

	args = vars(parser.parse_args())
	charger = args["charger"]
	maf = args["maf"]
	outBasename = args["outBasename"]
	outputDirectory = args["outputDirectory"]

	if outputDirectory[-1] != '/':
		outputDirectory = outputDirectory + '/'

	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	return charger, maf, outBasename, outputDirectory


###############
## MAIN CODE ##
###############

def main():
	charger, maf, outBasename, outputDirectory = argument_parser()

	# Open input file
	try:
		tsv=open(charger)
	except IOError:
		print("Input CharGer file does not exist!")

	# Open output file
	outFilename = outputDirectory+outBasename+".charged2vcf.filtered.tsv"
	outFilename_AF = outputDirectory+outBasename+".charged2vcf.filtered.af"+str(maf)+".tsv"
	
	outFile=open(outFilename, "w")
	outFileAF=open(outFilename_AF, "w")

	# Parse and filter CharGer file
	charger_header = tsv.readline()
	outFile.write(charger_header)
	outFileAF.write(charger_header)

	header = charger_header.strip().split("\t")

	for line in tsv:
		info = line.strip().split("\t")
		# filter out synonymous variants
		if info[header.index("Variant_Classification")] != "synonymous_variant":
			outFile.write(line)
			if float(info[header.index("Allele_Frequency")]) <= float(maf):
				outFileAF.write(line)

		# filter out variants with CharGer score < 4:
		if int(info[header.index("CharGer_Score")]) >= 4:
			outFile.write(line)
			if float(info[header.index("Allele_Frequency")]) <= float(maf):
				outFileAF.write(line)

	outFile.close()
	outFileAF.close()

if __name__ == "__main__":
	main()



