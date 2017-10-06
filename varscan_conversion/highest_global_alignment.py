#!/usr/bin/env python
import sys

sys.setrecursionlimit(1000000000)

indel_penalty = -10

#create score matrix for BLOSUM62 from file
def generate_blosum62_matrix():
	blosum62_matrix = {}
	
	blosum62_file = open("BLOSUM62.txt","r")

	first_line = next(blosum62_file).strip().split()

	for line in blosum62_file:
		if line != first_line:
			line = line.strip().split()
			columns = {}
			row_letter = ""
			for index,column in enumerate(line):
				if index == 0:
					row_letter = column
				else:
					columns[first_line[index - 1]] = int(column)
			blosum62_matrix[row_letter] = columns

	return blosum62_matrix

#used for debugging - print out matrix currently being used
def print_matrix(matrix, matrix_name):
	print "---------------------------------------------------------------"
	print "PRINTING " + matrix_name
	if matrix_name == "BLOSUM62":
		header_line = "-\t"
		for name in matrix.keys():
			header_line += name + "\t"
		print header_line
	for row in matrix:
		row_to_print = ""
		if matrix_name == "BLOSUM62":
			row_to_print += str(row) + "\t"
			row = matrix[row]
			for element in row:
				row_to_print += str(row[element]) + "\t"
		else:
			for element in row:
				row_to_print += str(element) + "\t"
		print row_to_print
	print "---------------------------------------------------------------"

#create Needleman-Wunsch matrix for alignment algorithm
def get_empty_NW_matrix():
	global string_one,string_two, traceback_matrix
	NW_matrix = []

	string_one = "-" + string_one
	string_two = "-" + string_two

	for first_index,first_letter in enumerate(string_one):
		row = []
		traceback_row = []
		for second_index,second_letter in enumerate(string_two):
			if first_index == 0:
				row.append(second_index * indel_penalty)
				
				if second_index == 0:
					traceback_row.append("done")
				else:
					traceback_row.append("left")
			
			elif second_index == 0:
				row.append(first_index * indel_penalty)
			
				if first_index != 0:
					traceback_row.append("up")
			else:
				row.append("-")
				traceback_row.append("unknown")
		NW_matrix.append(row)
		traceback_matrix.append(traceback_row)

	return NW_matrix

#determine which direction current location should be aligned in
def get_cell_score(row, column):
	global string_one, string_two, NW_matrix, blosum62_matrix, traceback_matrix

	diagonal_score = NW_matrix[row - 1][column - 1] + blosum62_matrix[string_one[row]][string_two[column]]
	up_score = NW_matrix[row - 1][column] + indel_penalty
	left_score = NW_matrix[row][column - 1] + indel_penalty

	max_dir = "diag"
	max_score = diagonal_score

	if left_score > max_score:
		max_score = left_score
		max_dir = "left"

	if up_score > max_score:
		max_score = up_score
		max_dir = "up"

	traceback_matrix[row][column] = max_dir
	return max_score

#move through two sequences and align
def fill_NW_matrix():
	global NW_matrix
	number_rows = len(NW_matrix)
	number_columns = len(NW_matrix[0])
	
	current_row_number = 0

	while current_row_number < number_rows:
		current_column_number = 0
		while current_column_number < number_columns:
			if NW_matrix[current_row_number][current_column_number] == "-": #meaning it has not been initialized
				score = get_cell_score(current_row_number, current_column_number)
				NW_matrix[current_row_number][current_column_number] = score
			current_column_number += 1
		current_row_number += 1

	return NW_matrix[number_rows - 1][number_columns - 1]

#move back through alignment and print out final words
def trace_path(row, column, final_string_one, final_string_two):
	global string_one, string_two, traceback_matrix

	movement = traceback_matrix[row][column] 

	if movement == "done":
		return [final_string_one[::-1],final_string_two[::-1]]
	elif movement == "diag":#match/inversion
		final_string_one += string_one[row]
		final_string_two += string_two[column]
		row -= 1
		column -= 1
	elif movement == "up":#deletion
		final_string_one += string_one[row]
		final_string_two += "-"
		row -= 1
	elif movement == "left":#insertion/duplication
		final_string_one += "-"
		final_string_two += string_two[column]
		column -= 1
	
	return trace_path(row, column, final_string_one, final_string_two)

#main method
def align(reference, alternate):
	global string_one, string_two, blosum62_matrix, traceback_matrix, NW_matrix
	
	string_one = reference
	string_two = alternate

	blosum62_matrix = generate_blosum62_matrix()
	
	traceback_matrix = []
	NW_matrix = get_empty_NW_matrix()
	best_score = fill_NW_matrix()

	return trace_path(len(string_one) - 1,len(string_two) - 1, "", "")
