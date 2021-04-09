from parse_structure import parse_structure
from H_bonds import find_Hbonds
from adapt_hbond import adapt_hbond


def cluster_hbonds(hbond):
    """
    Associate the H_bonds to a secondary structure
    as defined by the Dictionary of Secondary Structure (DSSP Kabsch & Sanders, 1983)
    """
    # Storing the starts and ends of secondary structures
    # in a dictionary
    # (to be extended for other SS)
    secondary_structure = {
        "3-turn": [], # [1, 8, 19, 30]
        "4-turn": [],
        "5-turn": []
    }

    # Coding for at least 3,4,5- turns

    # iterating over atoms indexes i and their corresponding neighbors j
    for i, j in enumerate(hbond):
        # i is not involved in a bond
        if (j != -1):

            # calculating the gap between neigbors
            n = abs(i - j)

            # For n-turns, n is in {3, 4, 5}
            if (n >= 3 and n <= 5):
                name = str(n) + "-" + "turn"

                # starting turn in i
                secondary_structure[name].append(i)

                n_new = n
                # are n_new and n of the same sign???
                while (n_new == n and j != -1):
                    i = j
                    j = hbond[i]
                    # comparing the gap with the neighbor of neighbor
                    n_new = abs(i - j)

                # ending turn in i + n
                secondary_structure[name].append(j)

    return secondary_structure


def test_clustering():

    alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens = parse_structure(
        "glut1.pdb")

    hbond = find_Hbonds(
        alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens)

    print("number of neighbors", len(adapt_hbond(hbond)))

    secondary_structure = cluster_hbonds(hbond)
    print("number of alpha helices", int(len(secondary_structure["4-turn"])/2))
    print("number of 3-helices", int(len(secondary_structure["3-turn"])/2))


test_clustering()
