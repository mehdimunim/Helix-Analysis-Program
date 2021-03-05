##################################################################################
# AUTHOR: Mehdi MUNIM
# Subroutines contains 80% of simulaid.f
# This program enumerates all of them and saves each subroutine in a separate file
# Assuming that subroutines are not nested
# Otherwise will return only the main subroutines
# The file has to be in the same directory as simulaid.f
###################################################################################

def parse_simulaid_f_by_subroutines():
	""" Create "subroutine XXX.f" for each simulaid's subroutine"""
	count_subroutines = 0
	with open("simulaid.f") as sf:
		# get each line of simulaid.f
		for line in sf:
			# get the beginning of a subroutine
			# startswith() avoids taking comment lines with 'subroutine'
			if (line.strip().startswith("subroutine")):
				count_subroutines +=1
				# save the "... subroutine ...." first line of the subroutine
				first_line = line
				title = first_line.strip().split('(')[0]
				print(title)
				# create a new FORTRAN subroutine file
				with open(title + ".f", "w+") as subroutine:
					# write the first line
					subroutine.write(first_line)
					# get the second line
					other_line = sf.readline()
					# copy each line of the subroutine's body in this new file
					# copy first the second line
					subroutine.write(other_line)
					# stop when "end" line is reach
					while other_line.strip() != "end":
						# copy the body of the subrouting in subroutine file
						other_line = sf.readline()
						subroutine.write(other_line)
	print("\nNumber of subroutines:",count_subroutines)


def main():
	parse_simulaid_f_by_subroutines()
main()
exit(0)
