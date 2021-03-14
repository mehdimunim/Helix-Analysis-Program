################################
# AUTHOR: Mehdi Munim
# Visualizes a static call graph
################################

import os


class Call_tree:

    def __init__(self, name, call_number):
        self.name = name
        self.call_number = call_number
        self.children = []

    def get_subroutines_call(self):
        """
        Get all subroutines list from a file
        """
        if (self.name == "simulaid.f"):
            file = self.name
        else:
            file = "subroutine " + self.name + ".f"
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
        return res

    def __eq__(self, tree):
        return self.name == tree.name

    def fill(self):
        """

        Fills the call tree recursively


        """
        # get every calls in the file
        subroutines = self.get_subroutines_call()
        # check if the node is orphane
        if (len(subroutines) != 0):
            for sub in subroutines:
                tree = Call_tree(sub, 1)
                self.children.append(tree)
                tree.fill()


def pprint_tree(node, file=None, _prefix="", _last=True):
    print(_prefix, "`- " if _last else "|- ", node.name, sep="", file=file)
    _prefix += "   " if _last else "|  "
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        _last = i == (child_count - 1)
        pprint_tree(child, file, _prefix, _last)


def main():
    main_tree = Call_tree("helixcomp", 1)
    main_tree.fill()
    pprint_tree(main_tree)


main()
exit(0)
