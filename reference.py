#! /usr/bin/evn python

import sys
import os

class ReferenceLoader:
    #constructor
    def __init__(self, correct_file, sv_type):
        self.correct_svs = open(correct_file, 'r')#tab-delimited file with information on correct svs
        self.sv_type = sv_type #type of structural variant we want to include in our reference
        
        self.population_positions = {}#positions of strain names in header; strain name -> column number
        self.variants_by_type = {}#variant_type maps to each variant position and strain information

        self.number_correct = 0
        self.create_reference()

    #correct svs file had complicated descriptions for the types of variations
    #this method simplifies to DEL, INS, INV, etc.
    def simplify_sv_name(self, name):
	simplified_name = name.upper()
	if("_" in name):
		full_name = name.split("_")
		simplified_name = full_name[1]
	elif("DEL" in name):
		simplified_name = "DEL"
	elif("INS" in name):
		simplified_name = "INS"
	elif("DUP" in name):
		simplified_name = "DUP"
	elif("INV" in name):
		simplified_name = "INV"
	elif("CNV" in name or "COPY" in name):
		simplified_name = "CNV"
	return simplified_name

    #go through tab-delimited file containing correct structural variants
    def create_reference(self):
	for line in self.correct_svs:
		parts = line.split()#separating by tab into columns

		if "CHR" == parts[0]:#This line is the header line. 
			#Getting names and positions of species strains
			number_populations = len(parts) - 1
			while number_populations > 4:#Four is the number of columns previous to population information
				self.population_positions[parts[number_populations]]=number_populations
				number_populations = number_populations - 1
		else:	
		    chromosome = parts[0]

	            position_start = int(parts[1])
		    position_end = int(parts[2])
		    position = [position_start,position_end]

		    variant_length = parts[3]
		    variant_type = self.simplify_sv_name(parts[4])
                    
                    if (self.sv_type == "ALL") or (variant_type == self.sv_type):
			self.number_correct += 1
                            
                        if (variant_type not in self.variants_by_type):
			    self.variants_by_type[variant_type] = []#adding variant type to dictionary if never seen before

			populations = {}
			
			for population in self.population_positions:
			    populations[population]  = False
			    if parts[self.population_positions[population]] == "1":
				populations[population] = True

			variant_position_strains = [chromosome, position, populations, variant_length]
			self.variants_by_type[variant_type].append(variant_position_strains)
	
	self.correct_svs.close()
