import math
import numpy as np

import os
os.chdir("../dependencies")
import dist2
from dependencies import *

def dssp(c, n1, n, nslt, line, index, inamcol1, inamcol2, iresncol1, iresncol2, nneig, ineig, nbox, indices, ihbneig, ixc, ixo, ixn, ixa, dssplab, idistdssp, ch, cres, enghb, iparal, iantiparal, nss, itypss, ifss, ilss, ires0, nconfig, iwdssp, iwrose, iwhead, ifail, radtodeg, maxrepconf, maxng, nnlistlen, maxbox, listlen, maxss, maxrsd, maxrec):
    """
    DSSP algorithm:
    Define Secondary Structure of Proteins
    Identify hydrogene-bonds and thus the alpha-helix

    ***
    Parameters:
    c:
    n1:
    n:
    nslt:
    line:
    index:
    inamcol1:
    inamcol2:
    iresncol1:
    iresncol2:
    nneig:
    ineig:
    nbox:
    indices:
    ihbneig:
    ixc:
    ixo:
    ixn:
    ixa:
    dssplab:
    idistdssp:
    ch:
    cres:
    enghb:
    iparal:
    iantiparal:
    nss:
    itypss:
    ifss:
    ilss:
    ires0:
    nconfig:
    iwdssp:
    iwrose:
    iwhead:
    ifail:
    radtodeg:
    maxrepconf:
    maxng:
    nnlistlen:
    maxbox:
    listlen:
    maxss:
    maxrsd: max residue
    maxrec:

    Returns:
    res

    """
    ing = 0
    ils = 0
    lnam = inamcol2-inamcol1+1
    ifail = 0
    iok = 1
    nres = 0
    nresok = 0
    iccprevf = 0
    icccurrf = 0
    icfound = 0
    iofound = 0
    infound = 0
    ihfound = 0
    iafound = 0

    # Processing the PDB file
    for ia in range[n1, nslt]:
        atnam[: lnam] = line[index[ia]][inamcol1: inamcol2]
        if (lnam == 4):
            atnam[5: 8] = ':   '
        if (lnam > 4):
            #TODO: where is leftadjustn ?
            leftadjustn(atnam, atnam, lnam)
        if (atnam[1: 4] == 'C   ' or atnam[1: 4] == ' C   '):
            icfound = 1
            ixc[nres+1] = ia
        elif (atnam[1: 4] == 'O   ' or atnam[1: 4] == ' O   '):
            iofound = 1
            ixo[nres+1] = ia
        elif (atnam[1: 4] == 'N   ' or atnam[1: 4] == ' N   '):
            infound = 1
            ixn[nres+1] = ia
        elif (atnam[1: 4] == 'CA  ' or atnam[1: 4] == ' CA  '):
            iafound = 1
            ixa[nres+1] = ia
        elif (atnam[1: 4] == 'H   ' or atnam[1: 4] == ' H   ' or atnam[1: 4] == ' D  ' or atnam[1: 4] == 'D   ' or atnam[1: 4] == 'HN  ' or atnam[1: 4] == ' HN '):
            ihfound = 1
            ch[0][nres+1] = c[0][ia].copy()
        if (ia == nslt):
            iresn = -1
        else:
            read(line(index(ia+1))(iresncol1: iresncol2), *, ERR=999) iresn
            if (ia == n1):
                ireso = iresn
        if (iresn != ireso):
            # New residue
            nres += 1
            iok = 1
            if (ireso != 0):
                if (icfound*infound*iofound < 1):
                    if (nconfig <= maxrepconf):
                        print('Residue {} is missing N, C or O'.format(ireso))
                    iok = 0
                elif (ihfound == 0):
                    # Generate H coordinates from C(prev), CA and N
                    iok = 0
                    ch[0][nres] = 0
                    if (iccprevf*iafound == 1):
                        dnh = 0.0
                        dnc = 0.0
                        for k in range(3):
                            rx[k] = 2.0*c[k][ixn[nres]]-c[k][ixa[nres]]-ccprev[k]
                            dnh = dnh+rx[k]**2
                            dnc = dnc+(c[k][ixn[nres]]-ccprev[k])**2
                        if (dnc < 4.0):
                            for k in range(3):
                                ch[k][nres] = c[k][ixn[nres]] + \
                                    rx[k]/math.sqrt(dnh)
                            iok = 1
                        else:
                            if (nconfig <= maxrepconf):
                                print(
                                    "Chain break at residue {} no H generated".format(nres))
                    else:
                        if (nconfig <= maxrepconf):
                            print("Could not generate H for residue'", ireso)
                if (icfound == 1):
                    iccprevf = 1
                    if (nres > 1 and ixn[nres] > 0):
                        if (dist2(ccprev, c[0][ixn[nres]]) > 4.0):
                            ihbneig[nres-1] = -1
                            if (nconfig <= maxrepconf):
                                print("Chain break between residues {} and {}".format(nres - 1, nres))
                    ccprev = c[0][ixc[nres]].copy()
                else:
                    iccprevf = 0
                    ccprev = np.zeros(3)
                if (iok == 1):
                    nresok += 1
                    ihbneig[nres] = 0
                else:
                    ihbneig[nres] = -1
            ireso = iresn
            icfound = 0
            iofound = 0
            infound = 0
            ihfound = 0
            iafound = 0
    if (nresok == 0):
        raise ValueError("ERROR: No residues containing C, O and N were found")
        ifail = 1
        return res
    threshold = -0.5/(0.42*0.20*332.0)
    for ir in range(nres):
        enghb[ir] = threshold+10**(-5)
        dssplab[ir] = ' '
        if (ixc[ir] > 0):
            cres[0][ir] = c[0][ixc[ir]].copy()
        else:
            cres[0][ir] = np.zeros(3)
    rchb = 9.2
    #TODO: where is nnlistsim function??
    nnlistsim(1, nres, cres, nneig, ineig, indices, nbox, rchb, ifail, maxng, nnlistlen, maxbox, listlen, 0)
    if (ifail > 0):
        print("Link-cell routine failed")
        for ir in range(nres):
            if (ihbneig[ir] >= 0):
                for jr in range(ir+3, nres):
                    if (ihbneig(jr) >= 0):
                        if (dist2(c[0][ixc[ir]], c[0][ixc[ir]]) < rchb**2):
                            eij1 = 1.0/math.sqrt(dist2(c[0][ixc[ir]], c[0][ixc[ir]]))
                            eij2 = 1.0/math.sqrt(dist2(c[0][ixc[ir]], ch[0][jr]))
                            eij3 = -1.0/math.sqrt(dist2(c[0][ixc[ir]], ch[0][jr]))
                            eij4 = -1.0/math.sqrt(dist2(c[0][ixc[ir]], c[0][ixc[ir]]))

                            eij= eij1+eij2+eij3+eij4;

                            eji1 = 1.0/math.sqrt(dist2(c[0][ixc[ir]], c[0][ixn[ir]]))
                            eji2 = 1.0/math.sqrt(dist2(c[0][ixc[ir]], ch[0][jr]))
                            eji3 = -1.0/math.sqrt(dist2(c[0][ixc[ir]], ch[0][jr]))
                            eji4 = -1.0/math.sqrt(dist2(c[0][ixc[ir]], c[0][ixn[ir]]))

                            eji = eji1+eji2+eji3+eji4

                            if (eij < enghb[ir]):
                                if (enghb[ir] == 0.0):
                                    write(6, 1157) ir, jr, ihbneig(ir), enghb(ir), eij
                                ihbneig[ir] = jr
                                enghb[ir] = eij
                            if (eji < enghb[jr]):
                                if (enghb[jr] == 0.0):
                                    write(6, 1157) jr, ir, ihbneig(jr), enghb(ir), eji
                                ihbneig[jr] = ir
                                enghb[jr] = eji
    else:
        for ir in range(nres):
            if (ihbneig[ir] >= 0):
                for jjr in range(nneig[ir]):
                    jr = ineig(jr, ir)
                    if (ihbneig[jr] >= 0 and jr > ir+2):
                        if (dist2(c[0][ixc[ir]], c[0][ixc[jr]]) < rchb**2):
                            
                            eij1 = 1.0/math.sqrt(dist2(c[0][ixo[ir]], c[0][ixn[jr]]))
                            eij2 = 1.0/math.sqrt(dist2(c[0][ixc[ir]], ch[0][jr]))
                            eij3 = -1.0/math.math.sqrt(dist2(c[0][ixo[ir]], ch[0][jr]))
                            eij4 = -1.0/math.sqrt(dist2(c[0][ixc[ir]], c[0][ixn[jr]]))

                            eij= eij1+eij2+eij3+eij4

                            eji1 = 1.0/math.sqrt(dist2(c[0][ixo[jr]], c[0][ixn[ir]]))
                            eji2 = 1.0/math.sqrt(dist2(c[0][ixc[jr]], ch[0][ir]))
                            eji3 = -1.0/math.sqrt(dist2(c[0][ixo[jr]], ch[0][ir]))
                            eji4 = -1.0/math.sqrt(dist2(c[0][ixc[jr]], c[0][ixn[ir]]))

                            eji = eji1+eji2+eji3+eji4

                            if (eij < enghb[ir]):
                                if (enghb[ir] == 0.0):
                                    print(" eHB update {} to {} eold = {} old partner = {} enew = {}".format(ir, jr, ihbneig[ir], enghb[ir], eij))
                                ihbneig[ir]=jr
                                enghb[ir]=eij
                            # end if
                            if (eji < enghb[jr]):
                                if (enghb(jr) == 0.0):
                                    print(" eHB update {} to {} eold = {} old partner = {} enew = {}".format(jr, ir, ihbneig[jr], enghb[ir], eij))
                                ihbneig[jr]=ir
                                enghb[jr]=eji
    # A hydrogen bond exist from C=O(i) to NH(ihbneig(i))
    nss=0
    iparal=np.zeros(nres)
    iantiparal=np.zeros(nres)
    ir=1
    while (ir < nres):
        # Look for the next SS element
        while (ir < nres and ihbneig[ir] <= 0):
            ir += 1
        nhbinc=ihbneig[ir]-ir
        irf=ir
        nss0=nss
        if (nhbinc > 2 and nhbinc <= 5):
            # Possible helix start
            ihfound=1
            nstep=1
            ir += 1
            while (ihfound == 1):
                nhbinc1=ihbneig[ir] - ir
                nhbinc2=ihbneig[ir+1]-(ir+1)
                nhbinc3=ihbneig[ir+2]-(ir+2)
                if (nhbinc1 != nhbinc):
                    if (nhbinc2 != nhbinc):
                        if (nhbinc1 == nhbinc2):
                            ihfound=0
                        elif (nhbinc3 != nhbinc):
                            ihfound=0
                if (ihfound == 1):
                    nstep += 1
                    irprev=ir
                    ir += 1
            if ((ir-irf) > 1):
                nss += 1
                ifss[nss]=irf
                ilss[nss]=ir+nhbinc-1
                itypss[nss]=nhbinc+3
            elif (nhbinc == -3):
                # Possible Lambda helix
                ihfound=1
                nstep=1
                ir=ir+1
                while (ihfound == 1):
                    nhbinc1=ihbneig[ir] - ir
                    if (nhbinc1 != nhbinc):
                        ihfound=0
                    if (ihfound == 1):
                        nstep += 1
                        irprev=ir
                        ir += 1
                if ((ir-irf) > 1):
                    nss += 1
                    ifss[nss]=irf
                    ilss[nss]=ir+nhbinc-1
                    itypss[nss]=9
                elif (ihbneig[ir] > 0):
                    # Possible sheet
                    isfound=1
                    isfirst=0
                    islast=0
                    npar=0
                    napar=0
                    ineigmax=0
                    ineigmin=nres
                    ifs=-1
                    while (isfound == 1):
                        jr=ihbneig[ir]
                        if (jr > 0):
                            if (ihbneig[jr] == ir):
                                # HB(i,j) = HB(j,i) => Antiparallel
                                iantiparal[ir]=1
                                ing=jr
                                ifs=ir-1
                                ils=ir+1
                        if (jr >= 2):
                            if (ihbneig[jr-2] > 0):
                                if (ihbneig[jr-2] == ir):
                                    # HB(j-1,i)=HB(i,j+1) => Parallel
                                    iparal[ir]=1
                                    ing=jr
                                    ifs=ir
                                    ils=ir+2
                        if (ir >= 1 and ir < nres):
                            jrm=ihbneig[ir-1]
                            if (jrm > 0):
                                if (ihbneig[jrm] == ir+1):
                                    # HB(i-1,j)=HB(j,i+1) => Parallel
                                iparal[ir]=1
                                ing=jrm
                                ifs=ir
                                ils=ir+2
                        else:
                            jrm=0
                        if (jrm > 2):
                            if (ihbneig[jrm-2] > 0):
                                if (ihbneig[jrm-2] == ir+1):
                                    # HB(i-1,j+1)=HB(j-1,i+1) => Antiparallel
                                    iantiparal[ir]=1
                                    ing=jrm
                                    ifs=ir
                                    ils=ir+2
                        if ((iparal[ir]+iantiparal[ir]) > 0):
                            if (isfirst == 0):
                                isfirst=ifs
                            if (ineigmin > ing):
                                ineigmin=ing
                            if (ineigmax < ing):
                                ineigmax=ing
                        elif (ir == 1):
                            isfound=0
                        elif (ir > 1):
                            if ((iparal[ir-1]+iantiparal[ir-1]) == 0):
                                isfound=0
                                islast=ils
                        npar=npar+iparal[ir]
                        napar=napar+iantiparal[ir]
                        ir+= 1
                        if ((islast-isfirst) > 2 and isfirst > 0):
                            nss += 1
                            ifss[nss]=isfirst
                            ilss[nss]=islast
                            if (npar*napar > 0):
                                itypss[nss]=3
                            elif (npar > 0):
                                itypss[nss]=1
                                if ((ineigmax-ineigmin) > (islast-isfirst)*2):
                                    itypss[nss]=4
                            elif (napar > 0):
                                itypss[nss]=2
                                if ((ineigmax-ineigmin) > (islast-isfirst)*2):
                                    itypss[nss]=5
                else:
                    ir += -1
                if (nss == maxss):
                    raise ValueError(
                        "ERROR; maximum number of secondary structure elements \({}\) has been reached redimension the program".format(maxss))
                    if (iwdssp > 0):
                        # write(iwdssp, 1002) maxss
                    ir=nres
                elif (nss == nss0):
                    # Neither helix nor sheet - just skip over
                else:
                    if (nss > 1):
                        if (ilss[nss-1] >= ifss[nss]):
                            ifss[nss]=ilss[nss-1]+1
                    else:
                        if (ifss[0] < 1):
                            ifss[0]=1
                    ir += 1
    for iss in range(nss):
        for ir in range(ifss[iss], ilss[iss]):
            dssplab[ir]=typc[itypss[iss]]
            idistdssp[itypss[iss]][ir]=idistdssp[itypss[iss]][ir]+1

    if (iwdssp > 0):
        if (nconfig <= 1):
            if (nss > 0):
                for i in range(nss):
                    print("SS# {} Residue index range: \({},{}\) Type: {}".format(
                        i, ires0 + ifss[i], ires0 + ilss[i], ssname[itypss[i]][1:lssname[itypss[i]]]))
                # write(iwdssp, 2005)(i, ires0+ifss(i), ires0+ilss(i), ssname(itypss(i))(1: lssname(itypss(i))), i= 1, nss)
            else:
                # write(iwdssp, 2006)
                print(" No secondary structure element was found")
        if (iwhead == 1):
            for i in range(9):
                print("line 1: residue number \(mod 10\)")
                print(
                    "lines 2-5: the digits of the residue number to which the {} residue number of line 1 is H-bonded \(if any\)".format(typc[i]))
                print(
                    "line 6: S for residues with bend angle \(CA[ir - 2] - CA[ir] - CA[ir+2]\) > 70 deg ")
                print("line 7: + or -, the signe of the ")
                print("CA[ir-1] - CA[ir] - CA[ir+1] - CA[ir +2] angle")
                print("line 8: secondary structure element type: ")
                # write(iwdssp, 2004)(typc(i), ssname(i)(1: lssname(i)), i = 1, 9)
        iresf=1
        while (iresf <= nres):
            iresl=min0[nres][iresf+49]
            for ic in range(50):
                charl1[ic]=' '
                charl2[ic]=' '
            for ir in range(max0[2][iresf], min0[iresl][nres-2]):
                if (ixa[ir-2]*ixa[ir]*ixa[ir+2] != 0):
                    ca1, ca2, cbend=angles(dist2(c[0][ixa[ir-2]], c[0][ixa[ir]]), dist2(c[0][ixa[ir+2]], c[0][ixa[ir]]), dist2(c[0][ixa[ir-2]], c[0][ixa[ir+2]]))
                    cosa=cbend
                    try:
                        bend=180.0 - dacoscheck(cosa)*radtodeg
                    except:
                        print("Error in dacoscheck angle=", cosa)
                    if (bend > 70.0):
                        charl1[ir-iresf+1]='S'
                else:
                    charl1[ir-iresf+1]='?'
            for ir in range(max0[1][iresf], min0[iresl][nres-2]):
                if (ixa[ir-1]*ixa[ir]*ixa[ir+1]*ixa[ir+2] != 0):
                    tors=dihangl(c, ixa[ir-1], ixa[ir], ixa[ir+1], ixa[ir+2])
                    if (tors >= 0.0):
                        charl2[ir-iresf+1]='+'
                    else:
                        charl2[ir-iresf+1]='-'
                else:
                    charl2[ir-iresf+1]='?'
            # write(iwdssp, 2000) iresf, iresl, (mod(i, 10), i=iresf, iresl)
            # write(iwdssp, 2001)(mod(ihbneig(i)/1000, 10), i=iresf, iresl)
            # write(iwdssp, 2001)(mod(ihbneig(i)/100, 10), i=iresf, iresl)
            # write(iwdssp, 2001)(mod(ihbneig(i)/10, 10), i=iresf, iresl)
            # write(iwdssp, 2001)(mod(ihbneig(i), 10), i=iresf, iresl)
            # write(iwdssp, 2002)(charl1(i-iresf+1), i=iresf, iresl)
            # write(iwdssp, 2002)(charl2(i-iresf+1), i=iresf, iresl)
            # write(iwdssp, 2002)(dssplab(i), i=iresf, iresl)
            iresf=iresl+1

    if (iwrose > 0):
        # Tentative alternative for turn detection (see Rose's 1977 paper)
        # write(iwrose, *) 'Data for turn detection with GW Rose method'
        rnprev=normplane(c[0][ixa[0]], c[0][ixa[2]], c[0][ixa[4]])
        rn1prev=normplane(c[0][ixa[1]], c[0][ixa[2]], c[0][ixa[3]])
        for ir in range(3, nres-2):
            r = radcirc(c[0][ixa[ir-2]], c[0][ixa[ir]], c[0][ixa[ir+2]])
            rn=normplane(c[0][ixa[ir-2]], c[0][ixa[ir]], c[0][ixa[ir+2]])
            rnn=np.dot(rn, rnprev)
            rn1=normplane(c[0][ixa[ir-1]], c[0][ixa[ir]], c[0][ixa[ir+1]])
            rnn1=np.dot(rn1, rn1prev)
            ca1, ca2, cbend=angles(dist2(c[0][ixa[ir-2]], c[0][ixa[ir]]), dist2(
                c[0][ixa[ir+2]], c[0][ixa[ir]]), dist2(c[0][ixa[ir-2]], c[0][ixa[ir+2]]))
            cosa=cbend
            try:
                bend=180.0-dacoscheck(cosa)*radtodeg
            except:
                print("Error in dacoscheck angle=", cosa)
            charl1[0]=' '
            if (bend > 70.0):
                charl1[0]='S'
            # write(iwrose, 2003) ir, charl1(1), dssplab(ir), r, rnn, rnn1
            rnprev=rn.copy()
            rn1prev=rn1.copy()
    raise ValueError(
        "ERROR: invalid residue number for atom {} : {}".format(ia, line))
    if (iwdssp > 0):
       # write(iwdssp, 1000) ia, line
    return res
