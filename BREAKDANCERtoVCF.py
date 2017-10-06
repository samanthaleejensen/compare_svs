#! /usr/bin/evn python

import sys
import os

BREAKDANCER_input = open(sys.argv[1], "r")
VCF_output = open(sys.argv[2], "w")

header_lines = "##fileformat=VCFv4.2\n\
##source=BREAKDANCER\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant detected\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF_output.write(header_lines)

for line in BREAKDANCER_input.readlines():
    if "#" not in line: # making sure it's not a header line
        line = line.split()
        chromosome = line[0]
        start = line[1]
        end = line[4]
        sv_type = line[6]
        length = abs(int(line[7]))

        vcf_line = chromosome + "\t" + start + "\t.\t.\t.\t.\tPASS\tSVTYPE=" + sv_type + ";SVLEN=" + str(length) + ";END=" + str(end) + "\n"

        VCF_output.write(vcf_line)

BREAKDANCER_input.close()
VCF_output.close()
