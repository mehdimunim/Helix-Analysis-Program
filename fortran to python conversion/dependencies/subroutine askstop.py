def askstop(idef):
    """
    ***
    Parameters:
    idef:
    """
    if (ipredict == 0):
        askyn('Do you want to continue', 23, 1, idef, icont, 20, 0)
        if (icont == 0):
            break
    else:
        print("Run continues as predictable input was requested")
    return
