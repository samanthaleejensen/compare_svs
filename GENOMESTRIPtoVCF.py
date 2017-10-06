#! /usr/bin/evn python

import sys
import os

GENOMESTRIP_input = open(sys.argv[1], "r")
VCF_output = open(sys.argv[2], "w")

header_lines = "##fileformat=VCFv4.2\n\
##source=GENOMESTRIP\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant detected\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF_output.write(header_lines)

for line in GENOMESTRIP_input.readlines():
    line = line.split()
    if "#" not in line[0]: # making sure it's not a header line
        chromosome = line[0]
        start = line[1]
        reference = line[3]
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
            elif "SVLEN=" in part:
                length = abs(int(part.replace('SVLEN=', '')))

        vcf_line = chromosome + "\t" + start + "\t.\t" + reference + "\t" + alternate + "\t.\t" + filter_value + "\tSVTYPE=" + sv_type + ";SVLEN=" + str(length) + ";END=" + str(end) + "\n"

        VCF_output.write(vcf_line)

GENOMESTRIP_input.close()
VCF_output.close()
