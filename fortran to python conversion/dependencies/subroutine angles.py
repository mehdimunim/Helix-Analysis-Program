import math


def angles(e1, e2, e3):
    """
    Calculate
    ***
    Parameters:
    e1:
    e2:
    e3:

    """
    ee1 = math.sqrt(e1)
    ee2 = math.sqrt(e2)
    ee3 = math.sqrt(e3)

    ca1 = (e2 + e3 - e1)/(2*ee2*ee3)
    ca2 = (e1 + e3 - e2)/(2*ee1*ee3)
    ca3 = (e1 + e2 - e3)/(2*ee1*ee2)
    return ca1, ca2, ca3
