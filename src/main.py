import sys
import os
import axis
import com
import common
import dssp as dssp_mod
import length
import parser
import tpr


def main():
    """
    Helix analysis program function to be launched directly


    """
    print("Helix Analysis Program")
    print("Greetings!")
    filename = sys.argv[1]

    print("Check arguments...")
    common.check_argument(sys.argv)

    print("Parsing file...")
    backbones, isTraj = parser.parse(filename, True)

    dssp = dssp_mod.DSSP(backbones)

    print("Saving helix assignements...")
    dssp.save_assignements()

    print("Turn angle per residue...")
    tpr.print_tpr(filename)

    if (isTrajectory):
        print("Length...")
        length.length_traj(filename)

    else:
        print("DSSP...")
        dssp.complete_print(filename)
        print("Inertia axes....")
        axis.inertia_axes(backbones)
        print("Length...")
        length.print_length(backbones)
        print("Center of mass...")
        com.COM(backbones)

    print("Done")


main()
