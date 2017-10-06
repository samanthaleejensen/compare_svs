#! /usr/bin/evn python

import sys
import os
from find_sv_type import find_sv_type

VARSCAN_input = open(sys.argv[1], "r")
VCF_output = open(sys.argv[2], "w")

header_lines = "##fileformat=VCFv4.2\n\
##source=VARSCAN\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant detected\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF_output.write(header_lines)

for line in VARSCAN_input.readlines():
    if "#" not in line: # making sure it's not a header line
        line = line.split()
        chromosome = line[0]
        start = line[1]
        
        # to calculate end, length, and sv_type, we need to know the alternate
        
        reference = line[3]
        alternates = line[4].split(",") # there may be more than one alternate allele in Varscan output
        
        for alternate in alternates:
            end = int(start) + len(alternate) - 1 #NOTE estimated - not reliable
            length = abs(len(alternate) - len(reference)) # again, estimated
            sv_type = find_sv_type(reference,alternate) # algorithm is not very robust - don't trust

            vcf_line = chromosome + "\t" + start + "\t.\t" + reference + "\t" + alternate + "\t.\tPASS\tSVTYPE=" + sv_type + ";SVLEN=" + str(length) + ";END=" + str(end) + "\n"

            VCF_output.write(vcf_line)

VARSCAN_input.close()
VCF_output.close()
