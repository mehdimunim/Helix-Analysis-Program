import pandas as pd
import numpy as np
import math
import os
from parse_structure import parse_structure
os.chdir(
    "C:\\Users\\Mehdi\\Documents\\GitHub\\Interdisciplinary-Project\\coding by ourselves")


def r(atom1, atom2):
    """
    Calculates distance between atom1 and atom2
    atom1, atom2 are 3-arrays of coordinates
    """
    return math.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)


def find_Hbonds(alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens):
    """
    For each residue, finds the nearest neigbor in term of H-Bonds energy.
    Returns the couple of residues engaged in H-Bonds
    """
    # total number of residues in the molecule
    number_residue = len(alpha_carbons)

    # Getting coordinates of the backbone chain
    backbone = {
        "oxygen": oxygens,
        "carbon": simple_carbons,
        "alpha_carbon": alpha_carbons,
        "nitrogen": nitrogens,
        "hydrogen": hydrogens
    }

    res = pd.DataFrame.from_dict(backbone)

    # Energy threshold to be considered as a hydrogen bond
    threshold = -0.5

    # Hbond_neighbor corresponds to the atoms engaged in an hydrogen bond
    # hbond[i] = j means that (i, j) are forming a bond
    hbond = [-1]*number_residue
    # the associated energies, initialized at threshold + 10^(-5)
    energy = [threshold + 10**(-5)]*number_residue

    # q1: partial charge in eV on the C of CO
    q1 = 0.42

    # q2: partial charge on the H of NH
    q2 = 0.20

    # dimensional factor
    f = 332

    # maximal tolerance radius between O(i) and N(j) for H-bonds
    # distance in angstoms (same as in the PDB)
    limit_distance = 5.2

    # Searching minimal energy H-bond in first direction
    for i in range(number_residue):
        # Criterium 0 to be an h-bond, having a distance of 3 residue

        for j in range(i + 3, number_residue):
            # Second criterium: residue-distance below limit
            distance = r(res["oxygen"][i], res["nitrogen"][j])
            if distance < limit_distance:

                # Bond energy in direction i -> j

                eijON = 1 / r(res["oxygen"][i], res["nitrogen"][j])
                eijCH = 1 / r(res["carbon"][i], res["hydrogen"][j])
                eijOH = 1 / r(res["oxygen"][i], res["hydrogen"][j])
                eijCN = 1 / r(res["carbon"][i], res["nitrogen"][j])

                eij = (q1 * q2 * f) * (eijON + eijCH - eijOH - eijCN)

                # Bond energy in the opposite direction j -> i

                ejiON = 1 / r(res["oxygen"][j], res["nitrogen"][i])
                ejiCH = 1 / r(res["carbon"][j], res["hydrogen"][i])
                ejiOH = 1 / r(res["oxygen"][j], res["hydrogen"][i])
                ejiCN = 1 / r(res["carbon"][j], res["nitrogen"][i])

                eji = (q1 * q2 * f) * (ejiON + ejiCH - ejiOH - ejiCN)

                # Second and third criterium: the two residues have the smallest energy and is below -0.5 kcal/mol
                if (eij < energy[i]):
                    energy[i] = eij
                    # i and j make up a minimal energy H-bond
                    hbond[i] = j

                if (eji < energy[j]):
                    energy[j] = eji
                    # i and j make up a minimal energy H-bond
                    hbond[j] = i
    return hbond
