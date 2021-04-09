def adapt_hbond(list):
    """
    Transforms the list "hbond" to a list of the couple of atoms in H-bonds
    Returns the reduced list
    """
    # if list[i] = -1, i is NOT in an H-bond
    value = - 1
    cleaned = []
    for index, item in enumerate(list):
        if (item != value):
            # +1 because residues index should start at 1
            cleaned.append((index + 1, item + 1))

    # eliminate duplicates
    # Not yet implemented

    return cleaned
