#! /usr/bin/evn python

import sys
import os
from find_sv_type import find_sv_type

#This class is what the FileParser method parse_next_line returns
#to keep all of the necessary information together in one place.
class Line:
    #Constructor for the object.
    def __init__(self, strain, chromosome, position_start, position_end, variant_type):
        self.strain = strain
        self.chromosome = chromosome
        self.position_start = position_start
        self.position_end = position_end
        self.variant_type = variant_type

    #ToString function for debugging.
    def printLine(self):
        print self.strain + " " + self.variant_type + ": " + self.chromosome + ", "  +  str(self.position_start) + "-" + str(self.position_end)

#This class contains all of the methods needed to parse a file and return the 
#tool_type, chromosome, position_start, position_end, and variant_type for 
#each line of a given input file
class FileParser:

    #This is the constructor for the class. It gathers all data about
    #the given file that we will need to parse the lines.
    def __init__(self, tool, file_name, sv_type):
        self.tool = tool
        self.file_name = file_name
        self.sv_type = sv_type
        self.output_file = self.get_output_file()
        self.current_line = 0
        self.strain_list = {}

    #This function uses command line input to determine the location
    #of the file.
    def get_output_file(self):
        output_file_name = ""
        while(output_file_name == ""):
            output_file_name = raw_input("Path to " + self.tool + " " + self.file_name + " file: ")
            if os.path.isfile(output_file_name) == False:
                print "Your entered file does not exist. Please try again."
                output_file_name = ""
        output_file = open(output_file_name, "r")
        return output_file.readlines()
    
    #This is the only function that will be called by other methods. It
    #returns the tool_type, chromosome, position_start, position_end, 
    #and variant_type for that line if it is the proper variant type.
    def parse_next_line(self):
        #make sure that we haven't reached the end of the file
        if self.current_line > len(self.output_file) - 1:
            return "EOF"

        next_line = self.output_file[self.current_line].split()
        self.current_line += 1
        
        #make sure it's not a header line
        if("#" not in next_line[0]) and ("CHROM" not in next_line[0].upper()) and (next_line[0] != "segStart"):
            strain = self.get_strain()
            chromosome = self.get_chromosome(next_line)
            position_start = self.get_start(next_line)
            position_end = self.get_end(next_line)
            sv_type = self.get_variant_type(next_line)

            #make sure it's the right kind of structural variant
            if(self.sv_type == sv_type) or (self.sv_type == "ALL"):
                return Line(strain, chromosome, position_start, position_end, sv_type)
        
        return self.parse_next_line()#it's a header line or the wrong sv type

    def get_strain(self):
        if(self.tool in ["DELLY","GENOME-STRIP","VARSCAN"]):#variant is from any of the included strains
            return "UNKNOWN"
        else:
            return self.file_name

    def get_chromosome(self, next_line):
        if (self.tool != "RDX"):
            return next_line[0]
        else:
            return next_line[7]

    def get_start(self, next_line):
        if (self.tool != "RDX"):
            return float(next_line[1])
        else:
            return float(next_line[8])

    def get_end(self, next_line):
        if (self.tool == "BREAKDANCER"):
            return float(next_line[4])
        elif (self.tool == "RDX"):
            return float(next_line[9])
        elif (self.tool == "CREST"):
            return float(next_line[5])
        else:#VCF format
            position_end = self.parse_extra_info(next_line[1], next_line[7])
            if (position_end == next_line[1]):#we didn't find the end position
                return self.get_bnd_end_position(next_line[1], next_line[4])
            return float(position_end)

    #the info VCF field holds the SV end position
    def parse_extra_info(self, position_start, excess_info):	
	excess_info = excess_info.split(';')
	for part in excess_info:
	    if "END=" in part and "CI" not in part:
		return int(part.replace('END=', ''))
	return position_start

    #if the variant doesn't have an SV end position field in a VCF
    #it may be in the alternate spot for some reason this seems to happen
    def get_bnd_end_position(self, position_start, possible_end_position):
	if ("<" not in possible_end_position) and ("<" not in possible_end_position):#it's a BND
        	end_position = possible_end_position.replace('[','')
		end_position = end_position.replace(']','')
		end_position = end_position.replace('N','')
		end_position = end_position.split(':')
		return int(end_position[1])
	return position_start
    
    def get_variant_type(self, next_line):
        if (self.tool == "VARSCAN") or (self.tool == "HYDRA-MULTI"):
        #these tools do not label their structural variants by type - to compare we need labeling capabilities
	    reference = next_line[3]
            alternate = next_line[4]
            alternates = alternate.split(",")#there could be more than one alternate allele per line
	    if len(alternates) > 1:
		output = find_sv_type(reference, alternates[0])
		solved = [alternates[0]]
		for possibility in alternates[1:]:
			if possibility not in solved:
				sv_type = find_sv_type(reference, possibility)
				output += "," + sv_type
				solved.append(possibility)
		return output
	    return find_sv_type(reference, alternate)
        elif (self.tool == "BREAKDANCER"):
            return next_line[6]
        elif (self.tool == "RDX"):
            state = next_line[2]
            if state == "1":
                return "DEL"
            elif state == "2":
                return "NONE"
            elif state == "3":
                return "DUP"
            else:
                return "UNKNOWN"
        elif (self.tool == "CREST"):
            return next_line[8]
        else:#VCF format
            return self.clean_vcf_variant_type(next_line[4])#will have <> or be gibberish

    #fixing vcf variant type to match notation in reference file
    def clean_vcf_variant_type(self, variant):
    	if '<' not in variant or '>' not in variant:
		return "BND"
	else:
		variant=variant.replace('<', '')
		variant = variant.replace('>', '')
		return variant
