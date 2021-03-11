def leftadjustline(line, ic1, ic2):
    """

    Adjust line on the left

    ***
    Parameters;
    line:
    ic1:
    ic2:

    Returns:
    line

    """
    icol = ic1
    icol = nextchar(line, icol, 132)
    nshift = icol - ic1
    if (nshift > 0):
        for i in range(icol, ic2):
            line[i-nshift:i-nshift+1] = line[i:i+1]
        for i in range(nshift):
            line[ic2-i+1:ic2-i+2] = " "
    return line
