#! /usr/bin/evn python

import sys

class VariantComparer:
    #constructor
    def __init__(self, allowable_margin, reference, output):
        self.allowable_margin = allowable_margin
        self.reference = reference
        self.output = output

        self.unique_variants = {}#will hold all unique calls made

        self.matches = 0
        self.near_matches = 0
        self.mislabeled = 0
        self.mismatches = 0
        self.variant_calls = 0


    #given two variants, determines if they are the same by comparing positions
    def variant_comparison(self, chromosome, position_start, position_end,  variant):
	whether_match = "NO MATCH"
	if (chromosome == variant[0]):
		start_position_difference = abs(variant[1][0] - position_start)
		end_position_difference = abs(variant[1][1] - position_end)
		if (start_position_difference <= self.allowable_margin) and (end_position_difference <= self.allowable_margin):
			whether_match = "MATCH"
		elif (start_position_difference <= 2*self.allowable_margin) and (end_position_difference <= 2*self.allowable_margin):
			whether_match = "NEAR MATCH"			
	return whether_match

    #determines if variant call is unique
    def is_unique(self, chromosome, position_start, position_end, variant_type):
	position = [position_start,position_end]#putting in the same format as variants_by_type so we can use variant_comparison method
	info = [chromosome, position]
	
	if variant_type not in self.unique_variants:
		self.unique_variants[variant_type] = []#adding variant type to dictionary if never seen before	
	else:
		for variant in self.unique_variants[variant_type]:
			if self.variant_comparison(chromosome, position_start, position_end, variant) == "MATCH":
				return False#this variant is already stored in the dictionary 	

	self.unique_variants[variant_type].append(info)
	self.variant_calls += 1	

	return True

    #for a particular variant type, is there any gold standard variant with a similar position?
    def get_matches_by_variant_type(self,test_var_type, variant_type, first_try, no_match, strain, chromosome, position_start, position_end):
	for variant in self.reference.variants_by_type[test_var_type]:	
		if(strain not in variant[2]):
			variant[2][strain] = "unknown"

		whether_match = self.variant_comparison(chromosome, position_start, position_end, variant)

		if whether_match != "NO MATCH":
			no_match = False
			if whether_match == "MATCH":		
				if first_try:#this is a true match
					self.matches +=  1

					if self.is_unique(chromosome, position_start, position_end, variant_type):	
						self.output.print_unique_match(variant, variant_type)
						self.output.print_variant(chromosome, position_start, position_end, variant_type, True)
				else:#the variant type was incorrectly labeled
					self.mislabeled = self.mislabeled + 1
					whether_match = "MISLABELED"
					
			elif whether_match == "NEAR MATCH":#within 2000 bp or incorrect strain	
				if first_try:
					self.near_matches += 1
				else:
					self.mislabeled += 1
					whether_match = "MISLABELED NEAR MATCH"

			self.output.print_variant_analysis(strain, whether_match, chromosome, position_start, position_end, variant_type, variant[0], variant[1][0], variant[1][1], test_var_type, variant[2][strain])
	return no_match

    #for a specific strain and variant, does there exist a matching gold-standard variant?
    def get_variant_results(self,strain, chromosome, position_start, position_end, variant_type):	
	global mismatches, variants_by_type
	whether_match = "NO MATCH"		
	no_match = True#since we can find multiple true and near matches, this variable holds whether we have found at least one already

	if variant_type in self.reference.variants_by_type:#checking for variant at position with correctly labeled variant type
		no_match = self.get_matches_by_variant_type(variant_type, variant_type, True, no_match, strain, chromosome, position_start, position_end)

	if no_match:#there are no true matches or near matches
		for var_type in self.reference.variants_by_type:#checking all SV types for position
			if variant_type != var_type: 
				no_match = self.get_matches_by_variant_type(var_type, variant_type, False, no_match, strain, chromosome, position_start, position_end)
		if no_match:#we found no matches, near matches, or mislabeled matches
			self.output.print_variant_analysis(strain, whether_match, chromosome, position_start, position_end, variant_type, "none", "none", "none", "none", "none")
			self.mismatches += 1
	
	if self.is_unique(chromosome, position_start, position_end, variant_type):#if it's a match it's already been checked
		self.output.print_variant(chromosome, position_start, position_end, variant_type, False)
	return

