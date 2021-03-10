import numpy as np
import math
import warnings


def multihelix(iw0, nhx, nhxres, radtodeg, c, icaahx, icbahx, icbreshx, maxat, maxhx, mxnhx):
    """

    Analyze multiple helix

    ***
    Parameters:
    iw0:
    nhx:
    nhxres:
    radtodeg:
    c:
    icaahx:
    icbahx:
    icbreshx:
    maxat:
    maxhx:
    mxnhx:

    Returns:
    res

    """
    MAXFRAMES = 50000
    MAXCOPY = 600
    MAXNHX = 12
    MAXNHX2 = MAXNHX2 * (MAXNHX2 - 1)
    ixres = nhx*nhxres+1
    print("MULTIHELIX nhxres =", nhxres)
    nframes = max0[0][nframe]
    nhx2tot = (nhx*(nhx - 1))/2
    i0line = 0
    line[1:80] = " "
    for ihx in range(nhxres):
        incr1 = (ihx-1)*nhxres
        ax1[0] = res[0][nframes][incr1+13]
        ax1[1] = res[1][nframes][incr1+13]
        ax1[2] = res[0][nframes][incr1+14]
        cent1[0] = res[1][nframes][incr1+14]
        cent1[1] = res[0][nframes][incr1+15]
        cent1[2] = res[1][nframes][incr1+15]
        halflen1 = res[0][nframes][incr1+7]/2.0
        for jhx in range(ihx+1, nhx):
            incr2 = (jhx-1)*nhxres
            ax2[0] = res[0][nframes][incr2+13]
            ax2[1] = res[1][nframes][incr2+13]
            ax2[2] = res[0][nframes][incr2+14]
            cent2[0] = res[1][nframes][incr2+14]
            cent2[1] = res[0][nframes][incr2+15]
            cent2[2] = res[1][nframes][incr2+15]
            halflen2 = res[0][nframes][incr2+7]/2.0
            ax1_ax2 = np.dot(ax1, ax2)
            try:
                ang = dacoscheck(dargin)*radtodeg
            except:
                print("Error in dacoscheck dargin = ", dargin)
            dist_cc = math.sqrt(dist2(cent1, cent2))
            c12 = cent1 - cent2
            c12_1 = np.dot(c12, ax1)
            c12_2 = np.dot(c12, ax2)
            parfact = 1.0
            if (abs(ax1_ax2) < 0.01):
                # Axes perpendicular
                a1 = -c12_1
                a2 = c12_2
                print(" Helices {} and {} are perpendicular", ihx, jhx)
            elif (abs(ax1_ax2) > 0.99):
                if (abs(ax1_ax2) > 0.999):
                    # Axes parallel
                    a2 = 0.0
                    ct1, a1 = project(cent2, cent1, ax1)
                    print("Helices {} and {} are parallel", ihx, jhx)
                    parfact = 0.1
                else:
                    for k in range(3):
                        dax1[k] = ax1[k]
                        dax2[k] = ax2[k]
                        dc12[k] = c12[k]
                    dax1_dax2 = np.dot(dax1, dax2)
                    dc12_1 = np.dot(dc12, dax1)
                    dc12_2 = np.dot(dc12, dax2)
                    a1 = (dc12_2-dc12_1/dax1_dax2) / \
                        (1.0/dax1_dax2-dax1_dax2)
                    a2 = (dc12_1-dc12_2/dax1_dax2) / \
                        (dax1_dax2-1.0/dax1_dax2)
            else:
                a1 = (c12_2-c12_1/ax1_ax2)/(1.0/ax1_ax2-ax1_ax2)
                a2 = (c12_1-c12_2/ax1_ax2)/(ax1_ax2-1.0/ax1_ax2)
            ct1 = cent1 + a1*ax1
            ct2 = cent2 + a2*ax2
            ct12 = ct1 - ct2
            ct12 = math.sqrt(np.dot(ct12, ct12))
            tst1 = np.dot(ct12, ax1)/c12norm
            tst2 = np.dot(ct12, ax2)/c12norm
            if ((tst1+tst2)*parfact > 0.01):
                if (tst1+tst2 < 0.1):
                    warnings.warns(
                        "Probable program Error: ihx = {} jhx = {} c12*m1 = {} c12*m2 = {} ".format(ihx, jhx, tst1, tst2))
                else:
                    raise Exception(
                        "Program Error: ihx = {} jhx = {} c12*m1 = {} c12*m2 = {} ".format(ihx, jhx, tst1, tst2))
            print(" I = {} JHX = {} Hlen1 = {} Hlen2 = {} a1 = {} a2 = {}".format(
                ihx, jhx, halflen1, halflen2, a1, a2))
            if (a1 < -halflen1 or a1 > halflen1 or a2 < -halflen2 or a2 > halflen2):
                # Make sure the closest approach points are inside the helix
                noproj = 0
                d2min = 99999.0
                s2 = cent2 - halflen2*ax2
                e2 = cent2 + halflen2*ax2
                s1 = cent1 - halflen1*ax1
                e1 = cent1 + halflen1*ax1
                exs1, as1 = project(s1, cent2, ax2)
                exe1, ae1 = project(e1, cent2, ax2)
                exs2, as2 = project(s2, cent1, ax1)
                exe2, ae2 = project(e2, cent1, ax1)
                if (d2min < 99999.0):
                    print(
                        " NOTE: closest HX {} - HX {} contact is outside", ihx, jhx)
                else:
                    if (noproj == 0):
                        print(
                            " NOTE: mismatch between  HX {} - HX {} ", ihx, jhx)
                    if (noproj == 1):
                        print(
                            "NOTE: no overlap between HX {} - HX {} ", ihx, jhx)
            angrot1 = 0
            if (icbahx[ihx] > 0):
                torsats[0][3] = c[0][icbahx[ihx]].copy()
                torsats[0][2], a = project(c[0][icbahx[ihx]], cent1, ax1)
                torsats[0][1] = torsats[0][2] + ax1
                s1, ac1 = project(cent2, cent1, ax1)
                ct12 = cent2 - s1
                ct12norm = math.sqrt(np.dot(ct12, ct12))
                for k in range(3):
                    ct12[k] /= ct12norm
                torsats[0][0] = torsats[0][1] + ct12
                angrot1 = dihangl(torsats, 1, 2, 3, 4, 1) * radtodeg
            angrot2 = 0
            if (icbahx[jhx] > 0):
                c[0][icbahx[jx]] = torsats[0][3].copy()
                torsats[0][2], a = project(c[0][icbahx[jhx]], cent2, ax2)
                torsats[0][1] = torsats[0][2] + ax2
                s2, ac1 = project(cent1, cent2, ax2)
                ct12norm = math.sqrt(np.dot(ct12, ct12))
                for k in range(3):
                    ct12[k] /= ct12norm
                torsats[0][0] = torsats[0][1] + ct12
                angrot2 = dihangl(torsats, 1, 2, 3, 4, 1)*radtodeg

            # Helix distances
            dist = math.sqrt(dist2(ct1, ct2))
            idist = math.floor(50*dist)
            idist_cc = math.floor(50*dist_cc)
            res[0][nframes][ixres] = 10000*idist + idist_cc

            # End-end distances
            idist_ss = math.floor(50*distss)
            idist_hh = math.floor(50*disthh)
            res[1][nframes][ixres] = 10000*idist_ss + idist_hh

            # Torsion angle over the closest approach "bond"
            torsats[0][1] = ct1
            torsats[0][2] = ct2
            torsats[0][0] = ct1 + halflen1*ax1
            torsats[0][3] = ct1 + halflen2*ax2
            dhang = dihangl(torsats, 1, 2, 3, 4, 1) * radtodeg
            iang = math.floor(10*ang)
            idhang = math.floor(10*dhang)
            res[0][nframes][nhx2tot + ixres] = 10000*iang+idhang

            # Torsion angle over the center-center "bond"
            torsats[0][1] = cent1.copy()
            torsats[0][2] = cent2.copy()
            torsats[0][0] = cent1 + ax1*halflen1
            torsats[0][3] = cent2 + ax2*halflen2
            dhang_cc = dihangl(torsats, 1, 2, 3, 4, 1)*radtodeg
            iamin1 = icaahx[icbreshx[ihx]][ihx]
            iamin2 = icaahx[icbreshx[jhx]][jhx]
            res[1][nframes][nhx2tot + ixres] = dhang_cc
            idist_ar1 = math.floor(10*angrot1)
            idist_ar2 = math.floor(10*angrot2)
            res[0][nframes][2*nhx2tot + ixres] = 10000*idist_ar1+idist_ar2
            print(" HXs# {} - {} dist= {} cc_dist = {} dee1 = {} dee2 = {} {}".format(ihx,
                                                                                      jhx, dist, dist_cc, distss, disthh, ap_pa(isg2ij+2)))
            print(" HXs# {} - {} ang= {} ang= {} dhang= {} dhang_cc = {} CAs: {} {}".format(
                ihx, jhx, ang, dhang, dhang_cc, iamin1, iamin2))
            print(" HX {} rotation wrt HX {} = HX {} rotation wrt HX {} = {}".format(
                ihx, jhx, angrot1, jhx, ihx, angrot2))
            ixres += 1
    return res
