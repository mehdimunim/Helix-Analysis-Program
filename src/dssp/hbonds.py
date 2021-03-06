#
# Finds the H-bonds
#

from ..common import distance as r      

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
    nresidue = len(alpha_carbon)

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
            distance = r(oxygen[i], nitrogen[j])
            if distance < limit_distance:

                # Bond energy in direction i -> j

                eijON = 1 / r(oxygen[i], nitrogen[j])
                eijCH = 1 / r(carbon[i], hydrogen[j])
                eijOH = 1 / r(oxygen[i], hydrogen[j])
                eijCN = 1 / r(carbon[i], nitrogen[j])

                eij = (q1 * q2 * f) * (eijON + eijCH - eijOH - eijCN)

                # Bond energy in the opposite direction j -> i

                ejiON = 1 / r(oxygen[j], nitrogen[i])
                ejiCH = 1 / r(carbon[j], hydrogen[i])
                ejiOH = 1 / r(oxygen[j], hydrogen[i])
                ejiCN = 1 / r(carbon[j], nitrogen[i])

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
