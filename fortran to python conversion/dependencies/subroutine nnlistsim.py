def nnlistsim(nfirst, n, c, nneig, ineig, indices, nbox, rcut, ifail, maxng, nnlistlen, maxbox, maxrec, LEVTEST):
    """


    ***
    Parameters:

    nfirst:
    n:
    c:
    nneig:
    ineig:
    indices:
    nbox:
    rcut:
    ifail:
    maxng:
    nnlistlen:
    maxbox:
    maxrec:
    LEVTEST

    Returns:
    Probably maxnn

    """
    # Set up neigbour list using linked cells
    if (n > nnlistlen):
        print("Redimension the program with maxat > ", 10*n):
        ifail = 1
        return
    div = rcut + 0.1
    rcut2 = rcut**2
    ifail = 0
    if (n < nfirst):
        return
    extension(c, nneig, 0, nfirst, n, xyzmin, xyzmax, centinp, 0, 0, v)
    nneig = n.zeros(n)
    for k in range(3):
        nxyz[k] = (xyzmax[k] - xyzmin[k])/div + 1
    nx = nxyz[0]
    nxy = nxyz[0]*nxyz[1]
    if (LEVTEST > 0):
        print("NNLISTSIM div={} nx,nxy={} {} nxyz= {} xyzmin= {} xyzmax= {}".format(
            div, nx, nxy, nxyz, xyzmin, xyzmax))
    nbox = np.zeros(ngrid)
    nboxmax = 0
    for i in range(nfirst, n):
        for k in range(3):
            ix[k] = (c[k][i] - xyzmin[k])/div + 1
        # Save i into the box represented by indices ix(1 - 3)
        indexi = ix[0] + nx*(ix[1] - 1) + nxy*(ix[2] - 1)
        nbox[indexi] = nbox[indexi] + 1
        if (nboxmax < nbox[indexi]):
            nboxmax = nbox[indexi]
        if (nbox[indexi] <= maxbox):
            indices[nbox[indexi]][indexi] = i
        else:
            ifail = 1
        if (ifail > 0):
            print("Too many atoms in a box increase maxbox and MAXREC and recompile")
            return None

    # Loop on the boxes
    for i1 in range(nxyz[0]):
        for i2 in range(nxyz[1]):
            for i3 in range(nxyz[2]):
                indexi = i1 + nx*(i2 - 1) + nxy*(i2 - 1)
                ni = nbox[indexi]
                if (LEVTEST > 2):
                    print(" NNLISTSIM i1,2,3= {}{}{} indexi= {} ni= {}".format(
                        div, nx, nxy, nxyz, xyzmin, xyzmax))
                if (ni > 0):
                    for j1 in range(max0[0][i1 - 1], i1):
                        for j2 in range(max0[0][i2 - 1], min0[nxyz[1], i2 + i1 - j1]):
                            j3lim = i3
                            if (i1 > j1 or i2 > j2):
                                j3lim += 1
                            for j3 in range(max0[0][i3-1], min0[j3lim][nxyz[2]]):
                                indexj = j1 + nx*(j2 - 1) + nxy*(j3 - 1)
                                nj = nbox[indexj]
                                # ni, nj are the number of atoms in the box (i1, i2, i3) and (j1, j2, j3)
                                if (i1 == j1 and i2 == j2 and i3 == j3):
                                    ij = 0
                                else:
                                    ij = 1
                                if (LEVTEST > 2):
                                    print("NNLISTSIM j1,2,3= {} {} {} indexj= {} nj= {}".format(
                                        j1, j2, j3, indexj, nj))
                                for jj in range(nj):
                                    j = indices[jj][indexj]
                                    ii1 = 1
                                    if (ij == 0):
                                        ii1 = jj + 1
                                    for ii in range(ii1, indexi):
                                        i = indices[ii][indexi]
                                        r2 = dist2(c[0][i], c[0][j])
                                        if (r2 <= rcut2):
                                            # Bond found
                                            if (nneig[i] < maxng and nneig[j] < maxng):
                                                nneig[i] += 1
                                                ineig[nneig[i]][i] = j
                                                neig[j] += 1
                                                ineig[nneig[j]][j] = i
                                            else:
                                                print(
                                                    "Redimension the program for more neighbors in nnlistsim")
                                                ifail = 1
                                                return None


maxnn = 0
for i in range(nfirst, n):
    if (nneig[i] > maxnn):
        maxnn = nneig[i]
return maxnn
