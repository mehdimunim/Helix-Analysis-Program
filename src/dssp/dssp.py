from .hbonds import find_hbonds
from .patterns import find_patterns


class DSSP:
    """
    Class implementing DSSP
    Inspired both from the paper of Kabsh & Sander and from Simulaid
    Contains print and get methods
    """

    def __init__(self, backbone):

        self.backbone = backbone

        self.hbonds = find_hbonds(self.backbone)

        self.secondary_structures, self.irreg = find_patterns(self.hbonds)

    def get_structures(self):
        return self.secondary_structures

    def get_helices(self):
        """
        Get list of start and end residus of all helices regardless of the type

        ---
        Return:

        all_helices: list of tuple (start, end)

        """
        all_helices = []
        for n in [3, 4, 5]:
            name = str(n) + "-helices"
            n_helices = self.secondary_structures[name]
            i = 1
            while (i < len(n_helices)):
                start = n_helices[i - 1]
                end = n_helices[i]
                all_helices.append((start, end))
                i += 2

        return all_helices

    def get_ca(self, res_number=False):
        """
        Get alpha carbons for each helix one by one
        ---
        Options:
        res_number: if true, will return lists of tuples of (coordinates, residue_number)

        ---
        Return:
        type: list of lists
        Each sub-list corresponds to a helix and contains the coordinates of its alpha-carbons
        """
        list = []
        res_num = self.backbone["res_number_list"]
        alpha_carbons = self.backbone["alpha_carbon"]

        for type in [3, 4, 5]:  
            name = str(type) + "-helices"
            helices = self.secondary_structures[name]
            i = 1
            while i < len(helices):

                # start and end residues of the helix
                start = helices[i - 1]
                end = helices[i]

                # alpha carbons between start and end
                if (res_number):
                    cas = [(alpha_carbons[j], res)
                           for j, res in enumerate(res_num) if res >= start and res <= end]
                else:
                    cas = [alpha_carbons[j]
                           for j, res in enumerate(res_num) if res >= start and res <= end]

                list.append(cas)

                i += 2

        return list

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

        print("\n **** Irregularities ****\n")
        for tuple in self.irreg:
            print(tuple)

    def complete_print(self, filename):
        """
        Print dssp output with as much information as possible
        """
        molecule_name = filename.split("/")[1][:-4]

        print("Molecule: ", molecule_name)
        secondary_structures = self.secondary_structures

        print("\nNumber of residues: ", len(self.backbone["alpha_carbon"]))
        print("\nStructures: ")
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
        print("\n Irregularities:\n")
        for tuple in self.irreg:
            print(tuple)

    def save_assignements(self, filename):
        """
        Saving DSSP's helix assignements to dssp.txt in folder output
        """
        molecule_name = filename.split("/")[1][:-4]
        output = "output/dssp_" + molecule_name + ".txt"

        secondary_structures = self.secondary_structures

        with open(output, "w+") as file:

            file.write("Molecule: " + molecule_name)
            nres = len(
                self.backbone["alpha_carbon"])
            file.write("\nNumber of residues: " + str(nres))
            file.write("\nStructures: ")

            for type in [3, 4, 5]:

                name = str(type) + "-helices"
                file.write(" \n {}\n".format(name.upper()))
                str_starts = ""
                str_ends = ""

                for pos, val in enumerate(secondary_structures[name]):
                    str_val = "{:4d}".format(val)
                    if pos % 2 == 0:
                        str_starts += str_val
                    else:
                        str_ends += str_val

                file.write("\nstarts: {:s} \n".format(str_starts))
                file.write("\nends:  {:s} \n".format(str_ends))
                file.write("\n Number of {}: {}\n".format(
                    name, int(len(secondary_structures[name])/2)))

            file.write("\n Irregularities:\n")

            for tuple in self.irreg:

                file.write(str(tuple))
