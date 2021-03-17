################################
# AUTHOR: Mehdi Munim
# Visualizes a static call graph
################################

import os


class Call_tree:

    def __init__(self, name):
        self.name = name
        # dictionary with children and number of calls
        self.children = {}

    def __eq__(self, other):
        return self.name == str(other)

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return (self.name)

    def get_subroutines_call(self):
        """
        Get all subroutines list from a file
        """
        # folder with all subroutines
        os.chdir(
            "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\subroutines\\all")
        # Exclude the main file simulaid.f
        if (self.name != "simulaid"):
            file = "subroutine " + self.name.strip()
        file += ".f"
        res = []
        if (os.path.isfile(file)):
            with open(file) as s:
                # get each line of simulaid.f
                for line in s:
                    # get the beginning of a subroutine
                    # startswith() avoids taking comment lines with 'subroutine'
                    if (line.strip().startswith("call")):
                        name = line.strip().split('(')[0][len("call"):]
                        res.append(name)
        else:
            print(file, "not found")
        return res

    def fill(self):
        """
        Fills the call tree recursively
        """
        # get every calls in the file
        subroutines = self.get_subroutines_call()
        # check if the node is orphane
        if (len(subroutines) != 0):
            for sub in subroutines:
                if sub in self.children:
                    self.children[sub] += 1
                else:
                    tree = Call_tree(sub)
                    tree.fill()
                    self.children.update({tree: 1})


def print_tree(node, prefix="", last=True, call_number=1):
    """
    Display the call graph on console
    Should also save it to a file
    Inspired from https://vallentin.dev/2016/11/29/pretty-print-tree
    """
    print(prefix, " \- " if last else "|- ", node.name,
          " (", call_number, ")", sep="")
    prefix += "   " if last else "|  "
    child_count = len(node.children)
    # put sort operations below
    for i, child in enumerate(node.children.keys()):
        last = i == (child_count - 1)
        call_number = node.children[child]
        print_tree(child, prefix, last, call_number)


def main():
    # Change the line below according to your desire
    main_tree = Call_tree("dssp")
    main_tree.fill()
    print_tree(main_tree)


main()
exit(0)
