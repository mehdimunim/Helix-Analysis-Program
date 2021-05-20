import sys
import os
import glob


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

    print("******** Helix Analysis Program ********\n")
    print("  Greetings!\n")

    # clear all files in output folder to avoid matplotlib error

    for file in glob.glob("output/*"):
        os.remove(file)
        print("Deleted " + str(file))

    print()
    print("  Checking arguments...\n")
    common.check_argument(sys.argv)

    filename = sys.argv[1]

    backbones, isTrajectory = parser.parse(filename, True)

    if (isTrajectory):

        print("   DSSP for the molecule in its initial position")

        dssp_init = dssp_mod.DSSP(backbones[0])

        dssp_init.complete_print(filename)

        print("   Saving helix assignements for the first frame...\n")
        dssp_init.save_assignements(filename)

        print("   Lengths of all helices during the trajectory...\n")
        length_mod.length_traj_corrected(filename)

        print("   Turn angle per residue for each helix...\n")
        tpr_mod.tpr_traj_corrected(filename)

    else:
        dssp = dssp_mod.DSSP(backbones)

        print("   Saving helix assignements...\n")
        dssp.save_assignements(filename)

        print("   DSSP...\n")
        dssp.complete_print(filename)

        print("   Length of helices...\n")
        length_mod.length(dssp.get_ca(), filename)

        print("   Turn angle per residue for helices...\n")
        tpr_mod.tpr(dssp.get_ca(), filename)

    print("Results are in output folder")
    print("Done")


main()
