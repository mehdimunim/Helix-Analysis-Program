import sys
import os


def main():
    """
    Helix analysis program function to be launched directly
    """
    sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..')))
    from src import length
    from src import tpr
    from src import axis
    from src import com
    from src import parser
    from src import common
    from src import dssp as dssp_mod

    print("Helix Analysis Program")
    print("Greetings!")

    filename = sys.argv[1]

    print("Check arguments...")
    common.check_argument(sys.argv)

    print("Parsing file {:s}...".format(filename))
    backbones, isTrajectory = parser.parse(filename, True)

    if (isTrajectory):
        print("Length...")
        length.length_traj(filename)

        print("Turn angle per residue...")
        tpr.print_tpr(filename)

    else:
        dssp = dssp_mod.DSSP(backbones)

        print("Saving helix assignements...")
        dssp.save_assignements(filename)

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
