#
# Finds the H-bond patterns as defined in Kabsch & Sander, 1983
#

def find_nturns(hbond):
    """
    Finds the beginnings of the n-turns with n in {3, 4, 5} with the list of H-bonds
    n-turns are H-bonds form CO(i) to NH(i+n)
    """
    n_res = len(hbond)
    nturn_starts = {
        "3-turn": [],
        # [1, 8] means there are two n-turns  (1, 1 + n) and (8, 8 +n)
        "4-turn": [],
        "5-turn": []
    }

    # iterating over the hbonds neigbors
    for bond in hbond:
        res1 = bond[0]
        res2 = bond[1]

        # Calculating the gap between H-bond neighbors
        n = res2 - res1

        if (n >= 3 and n <= 5):
            name = str(n) + "-turn"
            # Adding the beginning of the n-turn to the corresponding list
            nturn_starts[name].append(res1)

    return nturn_starts


def find_minimal_helices(nturn_starts):
    """
    Find helices on the basis of the n-turn beginnings lists
    Minimal helices are defined as consecutive n-turns
    """
    min_helices = {
        "3-min_helices": [],
        "4-min_helices": [],
        "5-min_helices": []
    }

    for n in [3, 4, 5]:

        name = str(n) + "-turn"
        list_nturns = nturn_starts[name]

        for i in range(len(list_nturns) - 1):
            if list_nturns[i+1] == list_nturns[i] + 1:
                helix_name = str(n) + "-min_helices"
                min_helices[helix_name].append(list_nturns[i])

    return min_helices


def assemble_minimal_helices(min_helices):
    """
    Gathers overlapping helices     
    """
    helices = {
        "3-helices": [],
        "4-helices": [],
        "5-helices": []
    }

    # gather the n-helices that start at less than n residues

    for n in [3, 4, 5]:
        input_name = str(n) + "-min_helices"
        output_name = str(n) + "-helices"

        list_min_helices = min_helices[input_name]
        n_res = len(list_min_helices)

        i = 0
        # avoiding special cases where there are 0 or 1 minimal helices
        if (n_res == 0):
            new_helix = []
        elif (n_res == 1):
            helices[output_name].append(list_min_helices[0])
            new_helix = []
        else:
            # starting new helix at position 0
            new_helix = [list_min_helices[0]]

        # iterate over consecutive indexes
        while (i < n_res - 1):
            # calculate the consecutive distance
            dist = list_min_helices[i + 1] - list_min_helices[i]

            # check if the helices overlap
            if (dist <= n):
                # if true add to current new helix
                new_helix.append(list_min_helices[i + 1])
            else:
                # else store the current helix
                helices[output_name].append(new_helix[0])
                helices[output_name].append(new_helix[-1] + n)

                # create a new helix at i + 1
                new_helix = [list_min_helices[i + 1]]

            # starting at next index
            i = i + 1

    return helices


def merge_helices(helices):
    """
    Merge helices that are consecutive, i.e less than 3 residues afar.
    Warns the user that there is an irregularities at the given precision 
    """
    irreg = []
    for n in [3, 4, 5]:
        name = str(n) + "-helices"
        n_helices = helices[name]
        i = 2
        while i < len(n_helices):
            # ending residue of the helix before
            end_last = n_helices[i-1]

            # beginning residue of the helix after
            start_new = n_helices[i]

            # number of residues between the helix before and after
            # note that gap > 0
            gap = start_new - end_last

            if gap <= 3:
                # helix irregularity between start_new and end_last
                irreg.append((end_last, start_new))

                # merge helices
                n_helices.remove(start_new)
                n_helices.remove(end_last)
            i += 2

    return helices, irreg


def find_patterns(hbond):
    """
    Associate the H_bonds to a secondary structure
    as defined by the Dictionary of Secondary Structure (DSSP Kabsch & Sanders, 1983)
    """

    n_turns = find_nturns(hbond)

    min_helices = find_minimal_helices(n_turns)

    helices = assemble_minimal_helices(min_helices)

    helices, irreg = merge_helices(helices)

    return helices, irreg
