import sys
import os


def main():
    """
    Helix analysis program function to be launched directly
    """
    sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..')))
    from src import length as length_mod
    from src import tpr as tpr_mod
    from src import axis as axis_mod
    from src import com as com_mod
    from src import parser
    from src import common
    from src import dssp as dssp_mod

    print("  Helix Analysis Program\n")
    print("  Greetings!\n")

    print("  Checking arguments...\n")
    common.check_argument(sys.argv)

    filename = sys.argv[1]

    backbones, isTrajectory = parser.parse(filename, True)

    if (isTrajectory):
        print("   Length...\n")
        length_mod.length_traj(filename)

        print("   Inertia axes...\n")
        axis_mod.inertia_traj(filename)

        print("   Turn angle per residue...\n")
        tpr_mod.tpr_traj(filename)

        print("   Center of mass...\n")
        com_mod.showGraphMassCenters(filename, True)

    else:
        dssp = dssp_mod.DSSP(backbones)

        print("   Saving helix assignements...\n")
        dssp.save_assignements(filename)

        print("   DSSP...\n")
        dssp.complete_print(filename)

        print("   Length...\n")
        length_mod.length(dssp.get_ca(), filename)

        print("   Inertia axes...\n")
        axis_mod.inertia(dssp.get_ca(), filename)

        print("   Turn angle per residue...\n")
        tpr_mod.tpr(dssp.get_ca(), filename)

        print("   Center of mass...\n")
        com_mod.showGraphMassCenters(filename, False)

    print("Done")


main()
