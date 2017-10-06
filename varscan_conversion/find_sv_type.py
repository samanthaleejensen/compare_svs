#! /usr/bin/evn python

import sys
import os

import highest_global_alignment

##This script can parse a line from the VCF and output the variant type field SVTYPE##
#########TODO fix duplication detection!!

#determines if added characters are arbitrary or duplications
##NOTE: if there is a duplication and another structural variant in the alternate 
##this will classify the whole variant as a duplication
def if_insertion_or_duplication(reference_alignment, alternate_alignment):
    if set(reference_alignment).discard("-") == set(alternate_alignment).discard("-"): #making sure they have only the same bases
        insertion_index = reference_alignment.find("-")
	
	current_index = insertion_index
	end = False	

	while end == False:
		if len(reference_alignment) == current_index + 1 or reference_alignment[current_index] != "-":
			end = True
		else:
			current_index += 1
		
	insertion = alternate_alignment[insertion_index:current_index]

	possible_duplications = [insertion]

	index = 1

	while index < len(insertion):
		if insertion[:index] in insertion[index:]:
			possible_duplications.append(insertion[:index])	
		index += 1
        
	for duplication in possible_duplications:
            if duplication != "":
		if reference_alignment[(insertion_index - len(duplication)):insertion_index] == duplication:
                        return "DUP"
		elif reference_alignment[current_index:(current_index + len(duplication))] == duplication:
                        return "DUP"
	
    return "INS"

#determines whether characters that differ are an inversion
def if_inversion_or_mismatch(reference_alignment, alternate_alignment):
	difference_indices = [i for i in range(len(reference_alignment)) if reference_alignment[i] != alternate_alignment[i]]
	
	reference_segment = reference_alignment[difference_indices[0]:difference_indices[-1]+1]
	alternate_segment = alternate_alignment[difference_indices[0]:difference_indices[-1]+1]

	if reference_segment == alternate_segment[::-1]:
		return "INV"
	else:
		return "BND"


def find_sv_type(reference, alternate):
	alignment = highest_global_alignment.align(reference, alternate)
	
	reference_alignment = alignment[0]
	alternate_alignment = alignment[1]

	if "-" in reference_alignment and "-" not in alternate_alignment:#possible insertion or duplication
		return if_insertion_or_duplication(reference_alignment, alternate_alignment)
	elif "-" not in reference_alignment and "-" in alternate_alignment:#deletion
		return "DEL"
	elif "-" not in reference_alignment and "-" not in alternate_alignment:#may be an inversion or simple mismatch
		return if_inversion_or_mismatch(reference_alignment, alternate_alignment)
	else:#super complicated and weird alignment
		return "BND"
		

