
def fit_angles(phi):
    """
    Fit theta such as:

    phi_i = theta * i + phi_0 (E)

    Solving the system:

        | SUM(E)
        | SUM(E*i for i)

    that can be written:

        | a11 * theta + a12 * phi_0 = b1
        | a21 * theta + a22 * phi_0 = b2

    ---
    Parameters:
    phi
    n

    ---
    Return:
    theta


    """
    n = len(phi)

    # coefficients
    a11 = (2*n + 1)*(n + 1)*n/6
    a21 = n*(n+1)/2
    a12 = a21
    a22 = n

    # Second member
    b1 = 0
    b2 = 0

    for i, phi in enumerate(phi):
        b1 += phi
        b2 += phi*i

    theta = (a22*b1 - a12 * b2)/(a22*a11 - a12*a21)

    return theta
