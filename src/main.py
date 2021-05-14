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

    print("  Helix Analysis Program\n")
    print("  Greetings!\n")

    filename = sys.argv[1]

    print("  Checking arguments...\n")
    common.check_argument(sys.argv)

    print("  Parsing file {:s}...\n".format(filename))

    backbones, isTrajectory = parser.parse(filename, True)

    if (isTrajectory):
        print("   Length...\n")
        length.trajectory_length(filename)

        print("   Inertia axes...\n")
        axis.inertia_axes(backbones)

        print("   Turn angle per residue...\n")
        tpr.print_tpr(filename)

    else:
        dssp = dssp_mod.DSSP(backbones)

        print("   Saving helix assignements...\n")
        dssp.save_assignements(filename)

        print("   DSSP...\n")
        dssp.complete_print(filename)

    print("Done")


main()
