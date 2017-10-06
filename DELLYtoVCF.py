#! /usr/bin/evn python

import sys
import os

DELLY_input = open(sys.argv[1], "r")
VCF_output = open(sys.argv[2], "w")

header_lines = "##fileformat=VCFv4.2\n\
##source=DELLY\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant detected\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF_output.write(header_lines)

for line in DELLY_input.readlines():
    if "#" not in line: # making sure it's not a header line
        line = line.split()
        chromosome = line[0]
        start = line[1]
        alternate = line[4]

        filter_value = line[6]

        #have to parse other fields for final elements
        sv_type = alternate.replace('<', '')
        sv_type = sv_type.replace('>', '')
       
        end = start
        length = "."

        info = line[7]
        info = info.split(';')
        for part in info:
            if "END=" in part and "CI" not in part:
                end = part.replace('END=', '')
            elif "CONSENSUS=" in part:
                alternate = part.replace('CONSENSUS=', '')
                length = len(alternate) # again, rough estimate <- not reliable

        vcf_line = chromosome + "\t" + start + "\t.\t.\t" + alternate + "\t.\t" + filter_value + "\tSVTYPE=" + sv_type + ";SVLEN=" + str(length) + ";END=" + str(end) + "\n"

        VCF_output.write(vcf_line)

DELLY_input.close()
VCF_output.close()
