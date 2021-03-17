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

    def get_subroutines_call(self):
        """
        Get all subroutines list from a file
        """
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

    def __eq__(self, other):
        return self.name == str(other)

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return (self.name)

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


def pprint_tree(node, file=None, _prefix="", _last=True):
    print(_prefix, "`- " if _last else "|- ", node.name,
          len(node.children.values()), sep="", file=file)
    _prefix += "   " if _last else "|  "
    child_count = len(node.children)
    for i, child in enumerate(node.children.keys()):
        _last = i == (child_count - 1)
        pprint_tree(child, file, _prefix, _last)


def print_tree(tree):
    return None


def main():
    main_tree = Call_tree("dssp")
    main_tree.fill()
    pprint_tree(main_tree)


main()
exit(0)
