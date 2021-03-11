import warnings


def extension(c, ih, nnh, nstart, n, cmin, cmax, c0, iprint, ivcheck, vol):
    """

    ***
    Parameters:
    c:
    ih:
    nnh:
    nstart:
    n:
    cmin:
    cmax:
    c0:
    iprint:
    ivcheck:
    vol

    """

    # Find an approximate center first
    for k in range(3):
        cmin[k] = 10**25
        cmax[k] = -10**25
    for k in range(3):
        if (nnh > 0):
            for ii in range(nstart, nnh):
                i = ih[ii]
                if (c[k][i] < cmin[k]):
                    cmin[k] = c[k][i]
                if (c[k][i] > cmax[k]):
                    cmax[k] = c[k][i]
        else:
            for i in range(nstart, n):
                if (c[k][i] < cmin[k]):
                    cmin[k] = c[k][i]
                if (c[k][i] > cmax[k]):
                    cmax[k] = c[k][i]
        c0[k] = (cmax[k] + cmin[k])/2
    vol = 1
    for k in range(3):
        vol *= (cmax[k] - cmin[k])
    if (iprint > 0):
        for k in range(3):
            print("Smallest middle and largest {} coordinate values = {} Volume of enclosing rectangle = {} A^3".format(
                xyz[k], cmin[k], c0[k], cmax[k], vol))
    if (ivcheck == 1:
        for k in range(3):
            imax=cmax[k]
            if (imax == 9999):
                warnings.warn(
                    "WARNING: input structure contains {} coordinate value(s) of 9999 indicating undefined values".format(xyz[k]))
        if (vol/n > 1000):
            warnings.warn(
                " WARNING: initial enclosing rectangle \(volume {} for {} atoms \) suggests corrupted input structure".format(vol, n))
            if (iprint == 0):
                for k in range(3):
                    print(" Smallest, middle and largest {} coordinate values= {} Volume of enclosing rectangle= {} A^3".format(
                        xyz[k], cmin[k], c0[k], cmax[k]))
            askstop(0)
    return
