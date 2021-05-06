#
# Finds the H-bonds
#

import math


def r(atom1, atom2):
    """
    Calculates distance between atom1 and atom2
    atom1, atom2 are 3-arrays of coordinates
    """
    return math.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)


def adapt_hbond(hbond, res_list):
    """
    Transforms hbond to a list of tuple (#residue, neighbor)
    """
    list = []
    for pos, val in enumerate(hbond):
        if val != - 1:
            list.append((res_list[pos], res_list[val]))
    return list


def find_hbonds(backbone):
    """
    For each residue, finds the nearest neigbor in term of H-Bonds energy.
    Returns the couple of residues engaged in H-Bonds
    """
    # Unpack backbone
    alpha_carbon = backbone["alpha_carbon"]
    carbon = backbone["carbon"]
    oxygen = backbone["oxygen"]
    nitrogen = backbone["nitrogen"]
    hydrogen = backbone["hydrogen"]
    res_number_list = backbone["res_number_list"]
    nresidue = len(res_number_list)

    # Energy threshold to be considered as a hydrogen bond
    threshold = -0.5

    # Hbond_neighbor corresponds to the atoms engaged in an hydrogen bond
    # hbonds[i] = j means that (i, j) are forming a bond
    hbonds = [-1]*nresidue

    # the associated energies, initialized at threshold + 10^(-5)
    energy = [threshold + 10**(-5)]*nresidue

    # q1: partial charge in eV on the C of CO
    q1 = 0.42

    # q2: partial charge on the H of NH
    q2 = 0.20

    # dimensional factor
    f = 332

    # maximal tolerance radius between N(i) and O(j) for H-bond
    # distance in angstoms (same as in the PDB)
    limit_distance = 5.2

    # Searching minimal energy H-bond in first direction
    for i in range(nresidue):
        # First criterium to be an h-bond, being at least 3 residue afar

        for j in range(i + 3, nresidue):
            # Second criterium: residue-distance below limit
            distance = r(backbone["oxygen"][i],
                         backbone["nitrogen"][j])
            if distance < limit_distance:

                # Bond energy in direction i -> j

                eijON = 1 / r(backbone["oxygen"][i], backbone["nitrogen"][j])
                eijCH = 1 / r(backbone["carbon"][i], backbone["hydrogen"][j])
                eijOH = 1 / r(backbone["oxygen"][i], backbone["hydrogen"][j])
                eijCN = 1 / r(backbone["carbon"][i], backbone["nitrogen"][j])

                eij = (q1 * q2 * f) * (eijON + eijCH - eijOH - eijCN)

                # Bond energy in the opposite direction j -> i

                ejiON = 1 / r(backbone["oxygen"][j], backbone["nitrogen"][i])
                ejiCH = 1 / r(backbone["carbon"][j], backbone["hydrogen"][i])
                ejiOH = 1 / r(backbone["oxygen"][j], backbone["hydrogen"][i])
                ejiCN = 1 / r(backbone["carbon"][j], backbone["nitrogen"][i])

                eji = (q1 * q2 * f) * (ejiON + ejiCH - ejiOH - ejiCN)

                # Second and third criterium: the two residues have the smallest energy and is below -0.5 kcal/mol
                if (eij < energy[i]):
                    energy[i] = eij
                    # i and j make up a minimal energy H-bond
                    hbonds[i] = j

                if (eji < energy[j]):
                    energy[j] = eji
                    # i and j make up a minimal energy H-bond
                    hbonds[j] = i

    return adapt_hbond(hbonds, res_number_list)


def test_hbonds():
    from parse import parse_structure
    backbone = parse_structure("glut1.pdb")

    hbonds = find_hbonds(backbone)
    print("#hbond", len(hbonds))
    print("example of hbond: ", hbonds[0])


# test_hbonds()
