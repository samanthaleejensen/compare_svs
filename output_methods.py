#! /usr/bin/evn python

import sys
import os

class OutputWriter:
    #constructor
    def __init__(self, output_name):
        #determine output names
        self.all_comparisons_name = output_name + ".all.txt"
        self.unique_name = output_name + ".unique.txt"
        self.sensitivity_name = output_name + ".matches.txt"
        self.log_name = output_name + ".log.txt"

        #open output files
        self.log_file = open(self.log_name, 'w')#output messages will be printed here
        self.all_comparisons_file = open(self.all_comparisons_name, 'w')#location positive predictive comparisons should be written
        self.unique_file = open(self.unique_name, 'w')#unique variant calls printed here
        self.sensitivity_file = open(self.sensitivity_name, 'w')#sensitivity comparisons will be written here

        #write headers
        self.all_comparisons_file.write("STRAIN\tMATCH TYPE\tCHR\tSTART\tEND\tTYPE\tCORRECT CHR\tCORRECT START\tCORRECT END\tCORRECT TYPE\tSTRAIN MATCHES\n")
        self.unique_file.write("CHR\tSTART\tEND\tTYPE\tMATCH\n")
        self.sensitivity_file.write("CHR\tSTART\tEND\tTYPE\tLENGTH\n")

        #initialize important values
        self.unique_matches = 0
        self.sv_length = 0

    #print to .match.txt file all unique true matches
    def print_unique_match(self, variant, variant_type):
	self.unique_matches = self.unique_matches + 1
	self.sv_length = self.sv_length + int(variant[3])
	self.sensitivity_file.write(variant[0] + "\t" + str(variant[1][0]) + "\t" + str(variant[1][1]) + "\t" + variant_type + "\t" + variant[3] + "\n")
	return

    #print all unique variant calls to .unique.txt file
    def print_variant(self, chromosome, position_start, position_end, variant_type, match_type):
	self.unique_file.write(chromosome + "\t" + str(position_start) + "\t" + str(position_end) + "\t" + variant_type + "\t" + str(match_type) + "\n")
	return

    #print output describing whether or not a match was found for each strain of a variant call
    def print_variant_analysis(self, strain, match_type, tool_chr, tool_start, tool_end, tool_type, correct_chr, correct_start, correct_end, correct_type, if_strain_matches):
	self.all_comparisons_file.write(strain + "\t" + match_type + "\t" + tool_chr + "\t" + str(tool_start) + "\t" + str(tool_end) + "\t" + tool_type + "\t" + correct_chr + "\t" + str(correct_start) + "\t" + str(correct_end) + "\t" + str(correct_type) + "\t" + str(if_strain_matches) + "\n")
	return

    #print desired message to screen and to log file
    def print_log(self, message):
	print message
	self.log_file.write(message + "\n")
	return

    def print_final_log(self, compare, reference):
        self.print_log("\n------TOOL ANALYSIS FINISHED------")
        self.print_log("Complete variant analysis by strain printed in " + self.all_comparisons_name + ".")
        self.print_log("Unique variant calls printed in " + self.unique_name  + ".")
        self.print_log("Unique matches printed in " + self.sensitivity_name + ".")
        self.print_log("Output log in " + self.log_name + ".")
        
        self.print_log("\nTOOL FOUND: " + str(compare.matches) + " matches, " + str(compare.near_matches) + " near matches, " + str(compare.mislabeled) + " incorrectly annotated svs, " + str(compare.mismatches) + " unmatched\n")

        if self.unique_matches != 0:
            self.print_log("Positive predictive value: " + str(self.unique_matches) + "/" + str(compare.variant_calls) + " (%.2f)" % (float(self.unique_matches)/float(compare.variant_calls)))
            self.print_log("Sensitivity: " + str(self.unique_matches) + "/" + str(reference.number_correct) + " (%.2f)" % (float(self.unique_matches)/float(reference.number_correct)))
            self.print_log("Average detected SV length: %.2f" % (float(self.sv_length)/self.unique_matches))
        
        self.log_file.close()
        self.all_comparisons_file.close()
        self.unique_file.close()
        self.sensitivity_file.close()

#OTHER METHODS

def print_usage():
    print "\n-----------------------------------------------------------------------\n\
ERROR: This script requires four  arguments - a file containing \n \
    the correct structural variants for the tested \n \
    region, the type of structural variant of interest, \n \
    the name of the SV discovery tool used, and a file\n \
    in which analysis will be outputted. \n \
\n\
USAGE: python compare_svs.py [CORRECT FILE] [SV TYPE] [TOOL NAME] [OUTPUT NAME]\n\
Ex: python compare_svs.py correct_svs.txt ALL lumpy lumpy_accuracy\n\
\n\
Type python compare_svs.py HELP for more information\n\
-----------------------------------------------------------------------\n"

def print_help(supported_tools):
	print "\n\t\t\t   COMPARE SVS HELP GUIDE\n\
-------------------------------------------------------------------------------\n\
\n\
This script compares output from various structural variant finders with known \n\
structural variants.It was designed for use with chromosome 19 of 8 common \n\
mouse strains and thus may not generalize well.\n\
\n\
-------------------------------------------------------------------------------\n\
1 - Usage\n\
-------------------------------------------------------------------------------\n\
This script requires four arguments: \n\
\n\
CORRECT FILE: a file of the format described below in section 3 which contains \n\
your desired standard reference, ideally a verified and complete list of \n\
structural variants in your region of interest. \n\
\n\
SV TYPE: an abbreviation for the variant type you'd like to compare (see \n\
section 4 for instructions). \n\
\n\
TOOL NAME: the name of the structural variant finder used to generate your \n\
data that you wish to compare to the standard (see section 5 for supported \n\
tools). \n\
\n\
OUTPUT NAME: your desired name for the output files that will be \n\
generated by this script. Keep in mind that any existing file of the same\n\
name will be replaced without warning. \n\
\n\
You should call the script as below:\n\
\n\
python compare_svs.py [CORRECT FILE] [SV TYPE] [TOOL NAME] [OUTPUT NAME]\n\
\n\
***NOTE*** \n\
In order to make this script friendly to those without computing skills \n\
it was designed with command line input of file locations required. \n\
Unfortunately this means that it is not able to be run as a supercomputing \n\
job. However, you will find that the comparisons and memory constraints \n\
of the script do not require great computing power and you should \n\
be able to run it from your personal computer if necessary. \n\
\n\
-------------------------------------------------------------------------------\n\
2 - Dependencies\n\
-------------------------------------------------------------------------------\n\
For this script to work, the following classes and methods must be found\n\
in the same directory as this main script, compare_svs.py: \n\
- parse_file_formats.py: contains scripts for parsing different file formats \n\
- output_methods.py: contains methods for printing to output files \n\
- reference.py: contains a class to hold all correct structural variants \n\
- variant_results.py: compare tool output to our reference structural variants \n\
- find_sv_type.py: used for SV tools that do not annotate structural variants \n\
- highest_global_alignment.py: a dependency of find_sv_type.py \n\
\n\
-------------------------------------------------------------------------------\n\
3 - Correct SVs File Format\n\
-------------------------------------------------------------------------------\n\
In order to create a dictionary of correct structural variants to compare \n\
results against, you will need to create a reference file in the same format \n\
as correct_svs.txt (should be found in this same directory).\n\
\n\
Each row in the file represents a different structural variant:\n\
COLUMN 1: chromosome \n\
COLUMN 2: start position (base pair) \n\
COLUMN 2: end position (base pair) \n\
COLUMN 3: variant length (base pairs) \n\
COLUMN 4: variant type* \n\
REMAINING COLUMNS: strains**\n\
\n\
*  Although correct_svs.txt has very specific variant types (ie Q6_del, H2_del) \n\
   this complexity is removed by the script that creates the reference dictionary \n\
   and is not necessary or even desired. Structural variant types will be reduced \n\
   to DEL, DUP, INV, INS, and CNV in most cases. \n\
\n\
** Each of these remaining columns represents whether or not the given structural \n\
   variant is present in the strain represented by the column header. A 0 here \n\
   indicates that the mutation is not found in that strain, while 1 means it is. \n\
   If you are attempting to use this script for just one species strain, simply \n\
   include only one column of whatever name you desire that is all ones. \n\
\n\
-------------------------------------------------------------------------------\n\
4 - Supported Variant Types\n\
-------------------------------------------------------------------------------\n\
For any of the variant types found in your reference file you can choose to \n\
compare your output just with that variant type. For example, with the \n\
correct_svs.txt file that this script was developed for, the fifth \n\
column contains complex variant type names that will be simplified to \n\
DEL, INV, and DUP. Those would be the options you could choose from. \n\
To compare all variant types, simply type ALL. \n\
\n\
-------------------------------------------------------------------------------\n\
5 - Supported Tools\n\
-------------------------------------------------------------------------------\n\
The following are structural variant finders whose output is currently \n\
interpretable by this script:"

	for tool in supported_tools:
		print "- " + tool

	print"\n\
If you have a tool not on this list that outputs VCF format results, \n\
you may be able to choose another VCF output tool (like LUMPY) \n\
and get accurate comparisons. \n\
\n\
To get another tool output format included email samleejensen@gmail.com. \n\
\n\
-------------------------------------------------------------------------------\n\
6 - Strains Tested\n\
-------------------------------------------------------------------------------\n\
The mouse strains in the dataset this script was developed for are as follows:\n\
- A_J\n\
- AKR_J\n\
- BALB_cJ\n\
- C3H_HeJ\n\
- C57BL_6NJ\n\
- CBA_J\n\
- DBA_2J\n\
- LP_J\n"
