#! /usr/bin/evn python

import sys
import os

CREST_input = open(sys.argv[1], "r")
VCF_output = open(sys.argv[2], "w")

header_lines = "##fileformat=VCFv4.2\n\
##source=CREST\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant detected\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF_output.write(header_lines)

for line in CREST_input.readlines():
    line = line.split()
    chromosome = line[18]
    start = line[19]
    end = line[22]
    sv_type = line[8]
    length = abs(int(line[20]) - int(line[17]) + 1) # estimate -> not reliable

    alternate = line[23] # from what I can tell this is the alternate -> not reliable, may also be reference (documentation is unclear)

    vcf_line = chromosome + "\t" + start + "\t.\t.\t" + alternate + "\t.\tPASS\tSVTYPE=" + sv_type + ";SVLEN=" + str(length) + ";END=" + str(end) + "\n"

    VCF_output.write(vcf_line)

CREST_input.close()
VCF_output.close()
