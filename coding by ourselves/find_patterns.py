from parse_structure import parse_structure
from H_bonds import find_Hbonds
from adapt_hbond import adapt_hbond


def find_nturns(hbond):
    """
    Finds the beginnings of the n-turns with n in {3, 4, 5} with the list of H-bonds
    n-turns are H-bonds form CO(i) to NH(i+n)
    """
    nturn_starts = {
        "3-turn": [],
        # [1, 8] means there are two n-turns  (1, 1 + n) and (8, 8 +n)
        "4-turn": [],
        "5-turn": []
    }

    # iterating over atoms indexes i and their corresponding neighbors j
    i = 0
    while i < n_res:

        # Calculating the gap between H-bond neighbors
        if (hbond[i] == - 1):
            n = - 1
        else:
            n = hbond[i] - i

        if (n >= 3 and n <= 5):
            name = str(n) + "-turn"
            # Adding the beginning of the n-turn to the corresponding list
            nturn_starts[name].append(i)
        i += 1

    return nturn_starts


def find_minimal_helices(nturn_starts):
    """
    Find helices on the basis of the n-turn beginnings lists
    Minimal helices are defined as consecutive n-turns
    """
    helices = {
        "3-helices": [],
        "4-helices": [],
        "5-helices": []
    }

    for n in [3, 4, 5]:
        name = str(n) + "-turn"
        list = sorted(nturn_starts[name])
        for i in range(len(list) - 1):
            if list[i+1] == list[i] + 1:
                helix_name = str(n) + "-helices"
                helices[helix_name].append(i)

    return helices


def assemble_minimal_helices(helices):
    """
    Gathers consecutive helices     
    """


def cluster_hbonds(hbond):
    """
    Associate the H_bonds to a secondary structure
    as defined by the Dictionary of Secondary Structure (DSSP Kabsch & Sanders, 1983)
    """
    n_turns = find_nturns(hbond)
    helices = find_helices(n_turns)


def test_clustering():

    alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens = parse_structure(
        "glut1.pdb")

    hbond = find_Hbonds(
        alpha_carbons, simple_carbons, oxygens, nitrogens, hydrogens)

    print("number of neighbors", len(adapt_hbond(hbond)))

    secondary_structure = cluster_hbonds(hbond)
    print("number of 3-helices", int(len(secondary_structure["3, 10-turn"])/2))
    print("number of alpha helices", int(len(secondary_structure["4-turn"])/2))
    print("number of pi helices", int(len(secondary_structure["5-turn"])/2))


test_clustering()
