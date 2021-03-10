import numpy as np
import math


def helixcomp(nslt, nres, calph0, dir0, perpvec0, start0, en0, cent0, rn0, camod, axfact, calph, dir, perpvec, start, end, cent, anglechange, anglechangeref, axtol, iw0, torsion, rotation, nrep, indexax, incrot, ireorienthx, idebughx, radtodeg, cc2, nreshx, icaahx, ihx, nhxres, maxhx, maxat):
    """
    Compute the rotation angle of a helix form the average angle between
    the corresponding perpendiculars to the helix axis from the Calphas
    ***
    Parameters:
    nslt:
    nres:
    calph0:
    dir0:
    perpvec0:
    start0:
    en0:
    cent0:
    rn0:
    camod:
    axfact:
    calph:
    dir:
    perpvec:
    start:
    end:
    cent:
    anglechange:
    anglechangeref:
    axtol:
    iw0:
    torsion:
    rotation:
    nrep:
    indexax:
    incrot:
    ireorienthx:
    idebughx:
    radtodeg:
    cc2:
    nreshx:
    icaahx:
    ihx:
    nhxres:
    maxhx:
    maxat:
    Returns

    """
    PI = math.pi
    incrhx = (ihx-1)*nhxres
    iverbort = 1
    if (nframe > 10 or nrep > 1):
        iverbort = 0
    print('HELIXC nframe,nrep,iverbort=', nframe, nrep, iverbort)
    if (nrep == 1 and nframe == 1 and idebughx > 0):
        print("initial helix cent:", cent0)
        for ir in range(nres):
            for k in range(3):
                print("Alpha C:", calph0[k][ir])
    it = 0
    # cc2 will be the possibly translated/overlaid frame
    res = helixaxis(cc2, nslt, 0, calph, dirw, startw, endw, cent, perpvec, camod, anglechangeref, circ, rn, axfact, axtol, rot, rms, helixlen, angles,
                    decidebend, nup, ndown, nrun, nnear, rcirc, turnperres, anglesn, 0, incrot, nrep, 1, nreshx, icaahx, ihx, nhxres, idebughx, radtodeg, pi, MAXHX)
    start = startw.copy()
    end = endw.copy()
    dir = dirw.copy()
    if (nrep == 1 and idebughx > 0):
        print("current helix cent:", cent)
        for ir in range(res):
            for k in range(3):
                print("Alpha C", calph[k][ir])
    # Shift the helix so that the start is at the origin
    for ir in range(nres):
        calph[0][ir] = calph[0][ir] - start
    end = end - start
    for k in range(3):
        shift[k] = start[k]
        start = 0
        startg = start.copy()
        endg = end.copy()
    if (nrep == 1):
        dirdotdir0 = np.dot(dir, dir0)
        try:
            rotation = dacoscheck(dirdotdir0)
        except:
            print("Error in dacoscheck rotation =", dirdotdir0)
        res[0][nframe][incrhx+5] = math.cos(rotation)
        res[1][nframe, incrhx+5] = sin(rotation)
        dev = cent - cent0
        res[0][nframe][incrhx+8] = math.sqrt(np.dot(dev, dev))
        res[1][nframe][incrhx+8] = dev[indexax[0]]
        res[0][nframe][incrhx+9] = dev[indexax[1]]
        res[1][nframe][incrhx+9] = dev[indexax[2]]
        devs = startw - start0
        res[0][nframe][incrhx+10] = devs[indexax[1]]
        res[1][nframe][incrhx+10] = devs[indexax[2]]
        deve = endw - en0
        res[0][nframe][incrhx+11] = deve[indexax[1]]
        res[1][nframe][incrhx+11] = deve[indexax[2]]
        if (idebughx > 0 and dirdotdir0 > 0.001):
            print("current helix cent", cent)
            for ir in range(nres):
                for k in range(3):
                    print("Alpha C:", calph[k][ir])
        if (ireorienthx == 1):
            # For better result, rotate dir onto dir0. The rotation axis is the
            # normal to the dir0-dir plane
            org = np.cross(dir0, dir)
            for k in range(3):
                corig[k][0] = 0.0
                ccurr[k][0] = 0.0
                corig[k][1] = dir0[k]
                ccurr[k][1] = dir[k]
                corig[k][2] = org[k]
                ccurr[k][2] = org[k]
            rot = ormat(ccurr, corig, 3, iverbort)
            for ir in range(nres):
                calph[0][ir] = dsmatvec(rot, calph(1, ir))
            dir = dir0.copy()
            if (idebughx > 0):
                print("current reoriented helix cent:", cent)
        for ir in range(nres):
            for k in range(3):
                print("Alpha C:", calph[k][ir])
    for ir in range(0, nres):
        perpvec[0][ir], org = calcperp(start, dir, calph[0][ir], it)
    nflip = 0
    rotav = 0
    rotav2 = 0
    changemin = 100.0
    changemax = 0.0
    nflipprev = nflip
    print("NRES-", nres)
    for ir in range(nres):
        angchange = angcomp(perpvec0[0][ir], dir0, perpvec[0][ir])
        print("ir={} dang = {} ".format(ir, angchange))
        for k in range(3):
            print("dir0 = {} perpvec = {} perpvec0 = {}".format(
                dir0[k], perpvec[k][ir], perpvec0[k][ir]))
        rotav += angchange
        rotav2 += angchange**2
        if (angchange < changemin):
            changemin = angchange
        if (angchange > changemax):
            changemax = angchange
        angchange[ir] = angchange

    torsion = rotav/nres
    sd = math.sqrt(abs(rotav2/nres-torsion**2))
    print("Nframe= {} avg = {} min = {} max = {}".format(
        nframe, torsion*180/PI, changemin*180/PI, changemax*180/PI))
    for ir in range(nres):
        print("anglechange[ir]*180/PI")
    if (changemax - changemin > PI/2):
        # Likely to have some sign flips
        for ir in range(nres):
            for k in range(3):
                xx[k] = -perpvec[k][ir]
            angchange = angcomp(perpvec0, dir0, xx)
            if (abs(angchange - torsion < abs(anglechange[ir] - torsion))):
                perpvec[0][ir] = xx
                nflip += 1
        if (nflip > 0):
            # Recalculate TPR
            turnperres + calcturnperres(nres, incrot,
                                        perpvec, dir, anglechangeref, 1, MAXHX)
            res[0][nframe][incrhx + 6] = math.cos(turnperres)
            res[1][nframe][incrhx + 6] = math.sin(turnperres)
        res[0][nframe][incrhx + 4] = math.cos(torsion)
        res[1][nframe][incrhx + 4] = math.sin(torsion)
        xx = np.cross(dir0, rn0)
        yy = np.cross(rn0, xx)
        for k in range(3):
            rot[0][k] = yy[k]
            rot[1][k] = xx[k]
            rot[2][k] = rn0[k]
        # Keep the normal from oscillating 180 degrees
        if (np.dot(rn, rn0) < 0):
            rn = dvmul(rn, -1)
            xx = dsmatvec(rot, rn)
            res[0][nframe][incrhx + 12] = xx[0]
            res[1][nframe][incrhx + 12] = xx[1]
            printhelix(iw0, startw, endw, cent, rms, helixlen, dirw, angles, decidebend,
                       nup, ndown, nrun, nnear, rcirc, turnperres, anglesn, ihx, radtodeg)
            try:
                rn_rn0 = radtodeg*dacoscheck(np.dot(rn, rn0))
            except:
                print("Error in dacoscheck angle=", np.dot(rn, rn0))
            print("Rotation = {} SD = {} Local tilt = {} N/Nr angle: {} X,Y,Z displacements = {} Absolute displacement = {}".format(
                radtodeg*torsion, radtodeg*sd, radtodeg*rotation, rn_rn0, dev, math.sqrt(np.dot(dev, dev))))
    if (nrep < 0):
        return radtodeg*torsion, radtodeg*sd, radtodeg*rotation, rn_rn0, dev, math.sqrt(np.dot(dev, dev))
    return radtodeg*torsion, radtodeg*sd, radtodeg*rotation, rn_rn0, dev, math.sqrt(np.dot(dev, dev))
