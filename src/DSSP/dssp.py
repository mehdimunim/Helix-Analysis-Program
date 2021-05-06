from parse import parse_structure
from hbonds import find_hbonds
from patterns import find_patterns


class DSSP:
    """
    Class implementing DSSP
    Inspired both from the paper of Kabsh & Sander and from Simulaid
    Contains print and get informations
    """

    def __init__(self, filename):

        self.backbone = parse_structure(filename)

        self.hbonds = find_hbonds(backbone)

        self.secondary_structures = find_patterns(hbond)

    def get_ca(self):
        """
        Get alpha carbons for each helix one by one
        """
        list_ca = []

        for type in [3, 4, 5]:
            name = str(type) + "-helices"
            helices = self.secondary_structures[name]
            for i in range(1, len(helices)):
                


        

        return list_ca

    def basic_print(self):
        secondary_structures = self.secondary_structures

        for type in [3, 4, 5]:

            name = str(type) + "-helices"
            print(" \n {}\n".format(name.upper()))

            str_starts = ""
            str_ends = ""
            for pos, val in enumerate(secondary_structures[name]):
                str_val = "{:4d}".format(val)
                if pos % 2 == 0:
                    str_starts += str_val
                else:
                    str_ends += str_val

            print("starts: ", str_starts)
            print("ends:   ", str_ends)

            print("Number of {}: {} ".format(
                name, int(len(secondary_structures[name])/2)))
