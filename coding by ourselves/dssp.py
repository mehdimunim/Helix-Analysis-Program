#  Criteria to defined secondary structures
# pattern recognition of hydrogen bonds
# Solvent exposure is the number of water molecules in contact with a residue
# RESULT:
# (1) Primary structure (SS bonds)
# (2) Secondary structure
# (3) Solvent exposure

# Two patterns to recognize
# (1) n-turns: H-bond between CO of residue i and the NH of residue i + n (n= 3, 4, 5)
# (2) bridges: H-bonds between residues not near each other in sequence
# alpha-helix => 4-turns

# Common imperfections: helical kinks, beta-bulges

# Structural overlaps solution: single state to each residue

# DEFINING AN H-BOND
# An ideal H bond has d= 2.9, theta = 0, E = - 3 kcal/mol
# MAX misalignement theta = 63
# E between H-bonding groups
# Partial charges on CO and NH
#                   1         1         1        1
# E = q1 * q2 (  ------- + ------- - ------- - -------  ) * f
#                 r(ON)     r(CH)     r(0H)     r(CN)
# q1 = 0.42e , q2 = 0.20e
# f= 332

# Minimal helix: two consecutive n-turns
# Longer helices are define as overlaps of minimal helices

# SECONDARY STRUCTURE IRREGULARITIES
# Proline kink

# BEND
# Region of high curvature

def parse_structure():
    """
    Parse the coordinates between: 
    Carbon atoms
    Oxygene atoms
    Nitrogen atoms
    Alpha carbons atoms
    Solvent atoms

    From it get the residues
    """
    # First step: removes solvent atoms

    # Second step: search for C, O, N and CA(alpha carbons)

    # Third step: get the different residues


def calculate_energies():
    """
    For each residue, finds the nearest neigbor in term of H-Bonds energy.
    Returns the energy values and the associated neighbour
    """
    # threshold to be considered as a hydrogene bond
    threshold = -0.5
    # q1: partial charge on the C of CO
    q1 = 0.42
    # q2: partial charge on the H of NH
    q2 = 0.20
    # dimensional factor
    f = 332
    # maximal tolerance radius for H-bonds
    rchb = 9.2
    for residue in atoms:
        for other_residue in range(residue + 3, number_of_atoms):
            if ( )
            # Comparing energies between CO and NH
            # in one side or in the other side
            eij1 = 1.0/r(c[0][ixo[ir]], c[0][ixn[jr]])
            eij2 = 1.0/r(c[0][ixc[ir]], ch[0][jr])
            eij3 = -1.0/r(c[0][ixo[ir]], ch[0][jr])
            eij4 = -1.0 / r(c[0][ixc[ir]], c[0][ixn[jr]])

            eij = (q1*q2*f)*(eij1+eij2+eij3+eij4)

            eji1 = 1.0/r(c[0][ixo[jr]], c[0][ixn[ir]])
            eji2 = 1.0/r(c[0][ixc[jr]], ch[0][ir])
            eji3 = -1.0/r(c[0][ixo[jr]], ch[0][ir])
            eji4 = -1.0/r(c[0][ixc[jr]], c[0][ixn[ir]])

            eji = (0.42*0.20*332.0)*(eji1+eji2+eji3+eji4)

            # Look for the smallest energy
            # so that for each residue, we only keep the minimal energy
            if (eij < enghb[ir]):
                # Neigborhood of the H-Bond
                neigborhood[ir] = jr
                energies[ir] = eij
            if (eji < energies[jr]):
                neigborhood[jr] = ir
                energies[jr] = eji


def find_helices():
    nss = 0
    iparal = np.zeros(nres)
    iantiparal = np.zeros(nres)
    ir = 1
    while (ir < nres):
        # Look for the next SS element
        while (ir < nres and ihbneig[ir] <= 0):
            ir += 1
        nhbinc = ihbneig[ir]-ir
        irf = ir
        nss0 = nss
        if (nhbinc > 2 and nhbinc <= 5):
            # Possible helix start
            ihfound = 1
            nstep = 1
            ir += 1
            while (ihfound == 1):
                nhbinc1 = ihbneig[ir] - ir
                nhbinc2 = ihbneig[ir+1]-(ir+1)
                nhbinc3 = ihbneig[ir+2]-(ir+2)
                if (nhbinc1 != nhbinc):
                    if (nhbinc2 != nhbinc):
                        if (nhbinc1 == nhbinc2):
                            ihfound = 0
                        elif (nhbinc3 != nhbinc):
                            ihfound = 0
                if (ihfound == 1):
                    nstep += 1
                    irprev = ir
                    ir += 1
            if ((ir-irf) > 1):
                nss += 1
                ifss[nss] = irf
                ilss[nss] = ir+nhbinc-1
                itypss[nss] = nhbinc+3
            elif (nhbinc == -3):
                # Possible Lambda helix
                ihfound = 1
                nstep = 1
                ir = ir+1
                while (ihfound == 1):
                    nhbinc1 = ihbneig[ir] - ir
                    if (nhbinc1 != nhbinc):
                        ihfound = 0
                    if (ihfound == 1):
                        nstep += 1
                        irprev = ir
                        ir += 1
                if ((ir-irf) > 1):
                    nss += 1
                    ifss[nss] = irf
                    ilss[nss] = ir+nhbinc-1
                    itypss[nss] = 9
                elif (ihbneig[ir] > 0):
                    # Possible sheet
                    isfound = 1
                    isfirst = 0
                    islast = 0
                    npar = 0
                    napar = 0
                    ineigmax = 0
                    ineigmin = nres
                    ifs = -1
                    while (isfound == 1):
                        jr = ihbneig[ir]
                        if (jr > 0):
                            if (ihbneig[jr] == ir):
                                # HB(i,j) = HB(j,i) => Antiparallel
                                iantiparal[ir] = 1
                                ing = jr
                                ifs = ir-1
                                ils = ir+1
                        if (jr >= 2):
                            if (ihbneig[jr-2] > 0):
                                if (ihbneig[jr-2] == ir):
                                    # HB(j-1,i)=HB(i,j+1) => Parallel
                                    iparal[ir] = 1
                                    ing = jr
                                    ifs = ir
                                    ils = ir+2
                        if (ir >= 1 and ir < nres):
                            jrm = ihbneig[ir-1]
                            if (jrm > 0):
                                if (ihbneig[jrm] == ir+1):
                                    # HB(i-1,j)=HB(j,i+1) => Parallel
                                iparal[ir] = 1
                                ing = jrm
                                ifs = ir
                                ils = ir+2
                        else:
                            jrm = 0
                        if (jrm > 2):
                            if (ihbneig[jrm-2] > 0):
                                if (ihbneig[jrm-2] == ir+1):
                                    # HB(i-1,j+1)=HB(j-1,i+1) => Antiparallel
                                    iantiparal[ir] = 1
                                    ing = jrm
                                    ifs = ir
                                    ils = ir+2
                        if ((iparal[ir]+iantiparal[ir]) > 0):
                            if (isfirst == 0):
                                isfirst = ifs
                            if (ineigmin > ing):
                                ineigmin = ing
                            if (ineigmax < ing):
                                ineigmax = ing
                        elif (ir == 1):
                            isfound = 0
                        elif (ir > 1):
                            if ((iparal[ir-1]+iantiparal[ir-1]) == 0):
                                isfound = 0
                                islast = ils
                        npar = npar+iparal[ir]
                        napar = napar+iantiparal[ir]
                        ir += 1
                        if ((islast-isfirst) > 2 and isfirst > 0):
                            nss += 1
                            ifss[nss] = isfirst
                            ilss[nss] = islast
                            if (npar*napar > 0):
                                itypss[nss] = 3
                            elif (npar > 0):
                                itypss[nss] = 1
                                if ((ineigmax-ineigmin) > (islast-isfirst)*2):
                                    itypss[nss] = 4
                            elif (napar > 0):
                                itypss[nss] = 2
                                if ((ineigmax-ineigmin) > (islast-isfirst)*2):
                                    itypss[nss] = 5


def check_irregularities():
     """Remove all proline residues from a pdb file
        Output : input pdb file without proline residues """
    
    with open("glut1.pdb", "r") as file_input:
    with open("glut1_removed.pdb", "w") as output: 
        for line in file_input:
               if line.startswith('ATOM'):
                    res_name = line[17:20].strip()
                    if res_name != 'PRO' :
                        output.write(line)
                        
                   
def print_output():


def dssp(atoms):
    parse_structure()
    calculate_energy()
    find_helices()
    check_irregularities()
    print_output()
