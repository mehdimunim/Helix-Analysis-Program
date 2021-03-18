######################################################
# AUTHOR: Mehdi Munim
# Convert a tree as defined in call_tree to DOT format
######################################################
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


def convert(tree, mypath):
    """
    Write tree structure in dot format 
    """
    os.chdir(mypath)
    # searching tree depth-first and writing the absolute path of the subroutines
    res = search_tree(tree)
    print(tree.name)
    with open(tree.name + ".dot", "w+") as dot:
        dot.write("graph " + tree.name + " { \n")
        for line in res:
            print(line)
            dot.write(line + "\n")
        dot.write("} \n")


def search_tree(tree):
    """
        Returns a list representation of tree in depth-first-search fashion
        each item is depicted with its absolute path from the mother node 
    """
    res = []

    def rec(tree, line):
        """
        Recursive function 
        """
        # if node is a leaf
        if len(tree.children) == 0:
            branch = "--".join(line)
            line = []
            res.append(branch)
        else:
            line.append(tree.name)
            for child in tree.children.keys():
                line.append(child.name.strip())
                rec(child, line)

    rec(tree, [])
    return res


def main():
    mypath = "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\subroutines\\call-tree-visualization"
    main_tree = Call_tree("dssp")
    main_tree.fill()
    convert(main_tree, mypath)


main()
exit(0)
