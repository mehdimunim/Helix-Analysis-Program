import numpy as np
import math


def kahn(co, nats, docircfit, iprint, message, MAXHX):
    """

    ref. Computers in Chemistry Vol 13, No 3, pg 191, 1989    

    Approach:    
        Construct a vector A from ca atom i to i-1. Construct B
    from i to i+1.  Find V1, the vector which bisects A and B.  V1 is
    perpendicular to the helix axis (for a perfect mathematical helix).
    Let i = i + 1 and repeat the procedure, finding V2.  As V1 and V2
    are both perp. to the helical axis, their cross product gives the
    helix direction.   Average over all possible tetrads of ca atoms,
    or use a 3-D fitting method to give the direction based on the
    points calculated to lie on the axis.    
    P1 and P2 are the vectors from the origin to CA 1 and CA 2.    
    The radius is a calculated as:    
        |dH|**2 - |p2-p1|**2
    r = --------------------
        2 * |(p1-p2) dot v2|    
    where d= (p2-p1) dot h    
    H1 and H2 are the position vectors
    H1 = P1 + r*V1
    H2 = P2 + r*V2    
    ***
    Parameters:
    ip: coordinates of the initial point
    fp: coordinates of the final point
    dir: the axis vector
    hsum: arrays of length MAXHX
    tmp, p1mp2, p2mp1, dmag, ddot: Arrays of length 3, 2*MAXHX    

    """
    for i in range(nats):
        for k in range(3):
            if (iprint > 3):
                print("Axis routine input:", co[k][i])
    hsum = np.zeros(3)
    hcount = 1
    for i in range(1, nats-2):
        # load p1 and p2
        p1 = co[0][i].copy()
        p2 = co[0][i+1].copy()
        # get vector v1
        a = p1 - co[0][i-1]
        b = p1 - co[0][i+1]
        a = dvnorm(a)
        b = dvnorm(b)
        v1 = a + b
        v1 = dvnorm(v1)
        # get vector v2
        a = p2 - co[0][i-1]
        b = p2 - co[0][i+1]
        a = dvnorm(a)
        b = dvnorm(b)
        v2 = a + b
        v2 = dvnorm(v2)
        # H= direction of axis
        h = np.cross(v1, v2)
        h = dvnorm(h)
        hsum += h
        # calculate radius
        p1mp2 = p1 - p2
        p2mp1 = p2 - p1
        # rise/residue
        d = np.dot(p2mp1, h)
        tmp = h*d
        # kahn
        r = (np.linalg.norm(tmp)**2 - np.linalg.norm(p2mp1)**2) / \
            (2*abs(np.dot(p2mp1, v2)))
        # calculate the points on the axis
        tmp = r*V1
        H1[0][hcount] = tmp + P1
        tmp = r*V2
        H1[0][hcount+1] = tmp + P2
        # hcount = hcount + 1
        hcount += 2
    dir = hsum.copy()
    dir = dvnorm(dir)
    dir, ip, fp, rms = parlsq(h1, hcount-1, docircfit, 0)
    if(docircfit):
        ip = circfit(co, nats, dir)
    else:
        # adjust ip to be next to the first alpha carbon, not the second
        tmp = co[0][0] - ip
        tmp = dvproj(dir, tmp)
        ip += tmp

    # seems to be a mistake because calculation
    #  of RMS and fp already done in parlsq (that uses RMScalc)
    # RMS, fp = RMScalc(co, nats, dir, ip)
    writeout_h(dir, ip, fp, rms, message, iprint)
    return rms, dir, ip, fp
