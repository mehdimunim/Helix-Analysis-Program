##############################################
# Author: Mehdi MUNIM
# Find every occurrences of words and 
# return the names of the corrresponding files 
##############################################

import os
import re

def find_files_with_words(word, path):
	""" find the files in path containing word
	return the files with the number of occurences found
	parameter word has to be in lowercase
	"""
	interesting_files = []
	# searching through file tree
	for root, dirs, files in os.walk(path):
		for subroutine in files:
			# searching for word in subroutine
			with open(subroutine) as sr:
				# putting the subroutine's content in lowercase
				code = sr.read().lower()
				# finding all occurences in the multiline file 
				res = re.findall(word,code, re.MULTILINE)
				# process when res is not null
				if res:
					# add the name of the file and the number of occurences found
					interesting_files.append( (subroutine,len(res)) )
	return interesting_files

def print_sorted(files_and_noc):
	""" print the subroutines by number of word occurences """
	for item in sorted(files_and_noc, key = lambda x : x[1], reverse = True):
		print(item[0], item[1])
	print("number: ", len(files_and_noc))


def main():
	word = "helix"
	# put the rights parameters in your case 
	files_and_noc = find_files_with_words(word, 'C:\\Users\\Mehdi\\EIDD\\Année 2\\Semestre 2\\Projet Transverse\\Simulaid\\subroutines\\subroutines')
	print_sorted(files_and_noc)

main()
exit(0)
