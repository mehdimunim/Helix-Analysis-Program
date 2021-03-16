################################
# AUTHOR: Mehdi Munim
# Visualizes a static call graph
################################

import os


class Call_tree:

    def __init__(self, name):
        # dictionary with the name of the subroutine and the number of occurences
        self.name = name
        self.name_and_occurences = {name: 1}
        self.children = []

    def get_subroutines_call(self):
        """
        Get all subroutines list from a file
        """
        # Exclude the main file simulaid.f
        if (self.name != "simulaid"):
            file = "subroutine " + self.name
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
                # adding one call else appending the new subroutine
                if sub in [child.name for child in self.children:
                    child.name_and_occurences[sub] += 1
                else:
                    tree = Call_tree(sub)
                    tree.fill()
                    self.children.append(tree)


def pprint_tree(node, file=None, _prefix="", _last=True):
    print(_prefix, "`- " if _last else "|- ", node.name,
          len(node.children), sep="", file=file)
    _prefix += "   " if _last else "|  "
    child_count = len(node.children)
    for i, child in enumerate(node.children):
        _last = i == (child_count - 1)
        pprint_tree(child, file, _prefix, _last)


def main():
    main_tree = Call_tree("dssp")
    main_tree.fill()
    pprint_tree(main_tree)


main()
exit(0)
