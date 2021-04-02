import pandas as pd
import numpy as np
import math


def r(atom1, atom2):
    """
    Calculates distance between atom1 and atom2
    atom1, atom2 are 3-arrays of coordinates
    """
    return math.sqrt(np.dot(atom1, atom2))


def find_Hbonds(coordinates, index_alpha_carbon, index_simple_carbon, index_oxygen, index_nitrogen, index_hydrogene):
    """
    For each residue, finds the nearest neigbor in term of H-Bonds energy.
    Returns the energy value and the associated neighbour
    """
    # Getting coordinates of the polypeptid skeleton
    res = {
        "oxygen": coordinates[index_oxygen],
        "carbon": coordinates[index_simple_carbon],
        "alpha_carbon": coordinates[index_alpha_carbon],
        "nitrogen": coordinates[index_nitrogen],
        "hydrogen": coordinates[index_hydrogene],

    }

    res = pd.DataFrame.from_dict(residues)
    # Dictionary for smallest neighbors index and corresponding energies
    hbond_neighbor = []
    energy = []

    # Energy threshold to be considered as a hydrogene bond
    threshold = -0.5

    # q1: partial charge on the C of CO
    q1 = 0.42

    # q2: partial charge on the H of NH
    q2 = 0.20

    # dimensional factor
    f = 332

    # maximal tolerance radius for H-bonds
    limit_distance = 9.2

    # Searching minimal energy H-bond in first direction
    for i in range(number_residue):
        # First criterium to be an h-bond, having a distance of 3 residue

        for j in range(i + 3, number_residue):
            # Second criterium: residue-distance below limit

            if (r(res["alpha_carbon", i]), res["alpha_carbon", j]) < limit_distance:

                eijON = 1 / r(res["oxygen", i], res["nitrogen", j])
                eijCH = 1 / r(res["carbon", i], res["hydrogen", j])
                eijOH = 1 / r(res["oxygen", i], res["hydrogen", j])
                eijCN = 1 / r(res["carbon", i], res["nitrogen", j])

                eij = (q1 * q2 * f) * (eijON + eijCH - eijOH - eijCN)

                # Third criterium: the two residues have the smallest energy
                if (eij < hbonds_neighbors[i]):
                    energy[i] = eij
                    # i and j make up a minimal energy H-bond
                    hbond_neighbor[i] = j

    # opposite direction
    for i in range(number_residue, 1, -1):
        # First criterium to be an h-bond, having a distance of 3 residue

        for j in range(i - 3, 1, -1):
            # Second criterium: residue-distance below limit

            if (r(res["carbon", i]), res["carbon", j]) < limit_distance:

                eijON = 1 / r(res["oxygen", i], res["nitrogen", j])
                eijCH = 1 / r(res["carbon", i], res["hydrogen", j])
                eijOH = 1 / r(res["oxygen", i], res["hydrogen", j])
                eijCN = 1 / r(res["carbon", i], res["nitrogen", j])

                eij = (q1 * q2 * f) * (eijON + eijCH - eijOH - eijCN)

                # Third criterium: the two residues have the smallest energy
                if (eij < hbonds_neighbors[i]):
                    # Neigborhood of the H-Bond
                    energy[i] = eij
                    hbond_neighbor[i] = j


def main():
    find_Hbonds()
