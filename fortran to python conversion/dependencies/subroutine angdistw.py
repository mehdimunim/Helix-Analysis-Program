def angdistw(c1, c2, c3, rhh):
    """
    Calculate the angle c2 - c1 - c3 and distance c1-c2; c1-c3
    Re
    """
    rroh1 = (c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2
    rroh2 = (c1[0]-c3[0])**2+(c1[1]-c3[1])**2+(c1[2]-c3[2])**2
    rrhh = (c2[0]-c3[0])**2+(c2[1]-c3[1])**2+(c2[2]-c3[2])**2
    droh1 = math.sqrt(rroh1)
    droh2 = math.sqrt(rroh2)
    rhh = math.sqrt(rrhh)
    cosa = (rroh1+rroh2-rrhh)/(2.d0*droh1*droh2)
    try:
        ahoh = dacoscheck(cosa)
    except:
        print("Problem in dacoscheck cosa = ", cosa)
    roh1 = droh1
    roh2 = droh2
    return roh1, roh2, rhh, ahoh
