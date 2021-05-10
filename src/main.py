import sys
import os
from parse import parse
from dssp import *
from print_tpr import print_tpr
from length_traj import *
from trajectory import inertia_axes
from print_length import print_length
from massCenter import COM

def main():
    """
    Helix analysis program function to be launched directly


    """
    print("Helix Analysis Program")
    print("Greetings!")
    filename = sys.argv[1]

    print("Check arguments...")
    check_argument(arguments)

    print("Parsing file...")
    backbones, isTraj = parse(filename, True)

    print("Saving helix assignements...")
    dssp.save_assignements()

    print("Turn angle per residue...")
    print_tpr(filename)

    if (isTrajectory):
        print("Length...")
        length_traj(filename)

    else:
        print("DSSP...")
        dssp.complete_print(filename)
        print("Inertia axes....")
        inertia_axes(backbones)
        print("Length...")
        print_length(backbones)
        print("Center of mass...")
        COM(backbones)

    print("Done")


main()
