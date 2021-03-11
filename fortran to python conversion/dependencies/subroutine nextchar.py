def nextchar(line, ifc, len):
    """

    Finds the next nonblank in line

    ***
    Parameters:
    line:
    ifc:
    len:

    Returns:
    ifc:

    """
    if (ifc > len - 1):
        return ifc
    ifc1 = ifc
    for i in range(ifc1, len-1):
        ifc = i
        if (line[i:i+1] != " " and line[i:i+1] != "\t"):
            return ifc
    return len
