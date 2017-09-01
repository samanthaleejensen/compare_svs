#! /usr/bin/evn python

import sys
import os

from parse_file_formats import *
from output_methods import *
from reference import ReferenceLoader
from variant_results import VariantComparer

###--------------------------------------------------COMPARE SVS--------------------------------------------------------###

##This script compares output from various structural variant finders with known structural variants.##
#It was designed for use with chromosome 19 of 8 common mouse strains and thus may not generalize well.#

#--------------------------------------------------------------------------------------------------------------------------
#DEPENDENCIES
#--------------------------------------------------------------------------------------------------------------------------

#   -parse_file_formats.py: contains scripts for parsing different file formats 
#   -output_methods.py: contains methods for printing to output files
#   -reference.py: contains a class that holds all correct structural variants
#   -variant_results.py: is used to compare tools to our reference structural variants
#   -find_sv_type.py: used for SV tools that do not annotate structural variants
#   -highest_global_alignment.py: a dependency of find_sv_type.py

#***All of the above files must be stored in the same directory as this script***

#---------------------------------------------------------------------------------------------------------------------------
#GLOBAL VARIABLES
#---------------------------------------------------------------------------------------------------------------------------

global output #an instance of the class OutputWriter from output_methods.py, does all log-writing and output printing
global reference #an instance of the class ReferenceLoader from reference.py, loads and holds reference structural variant file
global compare #an instance of the class VariantComparer from variant_results.py, compares a given structural variant to our reference variants
supported_tools = ["BREAKDANCER","LUMPY","RDX","CREST","DELLY","GENOME-STRIP","VARSCAN"]

#---------------------------------------------------------------------------------------------------------------------------
#HELPER FUNCTIONS
#---------------------------------------------------------------------------------------------------------------------------

#read through the given file of type FileParser
def parse_file(variant_file):
    found_variant = variant_file.parse_next_line()

    while(found_variant != "EOF"):#we haven't reached the end of the file yet
        compare.get_variant_results(found_variant.strain, found_variant.chromosome, found_variant.position_start, found_variant.position_end, found_variant.variant_type)
        found_variant = variant_file.parse_next_line()

#this method was added to parse the simulated heterozygote output since they had no true strains 
def get_strain_files(sv_finder, sv_type):
    global population_positions

    num_files = raw_input("Are your " + sv_finder + " files separated by strain? (Y/N) ")
    if num_files == "Y":
        for strain in population_positions:
            variant_file = FileParser(sv_finder, strain, sv_type)
            parse_file(variant_file)
    else:
        strains = raw_input("Enter the names of your strains separated by a slash (ex: A_J/AKR_J) ")
        variant_file = FileParser(sv_finder, strains, sv_type)
        parse_file(variant_file)

#----------------------------------------------------------------------------------------------------------------------------
#MAIN FUNCTION
#----------------------------------------------------------------------------------------------------------------------------

#check that there are the correct number of arguments (4)
if (len(sys.argv) == 2) and (sys.argv[1].upper() == "HELP"):
    print_help(supported_tools)
elif len(sys.argv) != 5:
    print_usage()
else:
    output = OutputWriter(sys.argv[4])

    #check that tool format is known and supported
    sv_finder_name = sys.argv[3].upper()#tool type output to compare is from
        
    if sv_finder_name not in supported_tools:
        output.print_log("ERROR: " + sv_finder_name + " not supported yet. \nSupported formats are: ")
        for tool in supported_tools:
            output.print_log("\t-" + tool)
        output.print_log("Contact Samantha Jensen at samleejensen@gmail.com to get new format included.")
	sys.exit()

    #check that structural variant type is present in reference file
    sv_type = sys.argv[2].upper()#type of structural variant detected that we want to compare
    reference = ReferenceLoader(sys.argv[1], sv_type)

    if len(reference.variants_by_type) == 0:
        print_log("ERROR: " + sv_type + " is not found in the given reference file.\n\
               Please choose another variant type to analyze.")
        sys.exit()

    output.print_log("------CORRECT VARIANTS LOADED SUCCESSFULLY------")
    output.print_log(str(reference.number_correct) + " correct variants of type " + sv_type)
    output.print_log("\nLoading tool output format...")
    
    compare = VariantComparer(1000, reference, output)

    if sv_finder_name == "DELLY":
        for variant_type in reference.variants_by_type:
            variant_file = FileParser(sv_finder_name, variant_type, variant_type)
            parse_file(variant_file)            
    elif sv_finder_name == "GENOME-STRIP" or sv_finder_name == "VARSCAN":
        variant_file = FileParser(sv_finder_name, "combined output", sv_type)
        parse_file(variant_file)
    else: #tools with separate strain files
        get_strain_files(sv_finder_name, sv_type)

    output.print_final_log(compare, reference)

###------------------------------------------------------------------------------------------------------------------------###

#NOTE to add another SV tool output to this script, see the file parse_file_formats.py. 
#Email samleejensen@gmail.com if you have questions.
