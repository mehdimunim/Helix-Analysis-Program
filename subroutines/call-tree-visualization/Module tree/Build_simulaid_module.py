#####################################################
# Transform the subroutines of simulaid.f into a module
# in order to implement fortrancallgraph
#  https://github.com/fortesg/fortrancallgraph
#####################################################
import os


def parse_simulaid_f_by_subroutines():
    """
    Creates "subroutine XXX.f" for each simulaid's subroutine
    and returns the list of all subroutines
    """
    subroutines = []
    with open("simulaid.f") as sf:
        # get each line of simulaid.f
        for line in sf:
            # get the beginning of a subroutine
            # startswith() avoids taking comment lines with 'subroutine'
            if (line.strip().startswith("subroutine")):
                # save the "... subroutine ...." first line of the subroutine
                first_line = line
                title = first_line.strip().split('(')[0]
                print(title)
                subroutines.append(title + ".f")
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
    return subroutines


def build_simulaid_module():
    subroutines = parse_simulaid_f_by_subroutines()
    with open("simulaid_module.f", "w+") as module:
        module.write("module simulaid\n")
        module.write("contains\n")
        for subroutine in subroutines:
            with open(subroutine, "r") as sub:
                for line in sub:
                    module.write("    " + line)
            os.remove(subroutine)
        module.write("\nend module simulaid")


def main():
    # remove the line below if necessary
    os.remove("simulaid_module.f")
    build_simulaid_module()


main()
exit(0)
