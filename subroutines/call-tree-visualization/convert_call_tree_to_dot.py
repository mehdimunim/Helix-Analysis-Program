######################################################
# AUTHOR: Mehdi Munim
# Convert a tree as defined in call_tree to DOT format
######################################################
import os
from copy import copy


class Call_tree:

    def __init__(self, name):
        self.name = name.strip()
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
        if (self.name == "simulaid"):
            file = "simulaid.f"
        else:
            file = "subroutine " + self.name.strip() + ".f"
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


def convert(tree, mypath):
    """
    Write tree structure in dot format
    """
    os.chdir(mypath)
    # searching tree depth-first and writing the absolute path of the subroutines
    res = search_tree(tree)
    with open(tree.name + ".dot", "w+") as dot:
        dot.write("digraph " + tree.name + " { \n")
        for line in res:
            dot.write(line + "\n")
        dot.write("} \n")


def search_tree(tree):
    """
        Returns a list representation of tree in depth-first-search fashion
        each item is depicted with its absolute path from the mother node
    """
    res = []

    def rec(node, line):
        """
        Recursive function that fills res
        """
        # if node is a leaf
        if len(node.children) == 0 and len(line) != 0:
            # cleaning the list of the elements in double
            start = []
            start.append(line[0])
            end = []
            end.append((line[-1][0], 0))
            line = line[1:-1]
            line = line[::2]
            line = start + line + end
            branch = "-->".join(
                [item[0] for item in line])
            res.append(branch)
        else:
            line.append((node.name, ""))
            for child in node.children:
                line.append((child.name, node.children[child]))
                rec(child, line.copy())
                # removes the name of the node that has been visited
                line.pop(-1)

    rec(tree, [])
    return res


def main():
    mypath = "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\subroutines\\call-tree-visualization"

    list_of_subroutines = ["analyze", "dssp",
                           "multihelix", "helixcomp", "checkforhelix", "parlsq", "kahn", "simulaid"]
    for sub in list_of_subroutines:
        main_tree = Call_tree(sub)
        main_tree.fill()
        convert(main_tree, mypath)


main()
exit(0)
