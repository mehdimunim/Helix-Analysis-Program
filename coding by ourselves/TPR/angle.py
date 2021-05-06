
def angle(v1, v2, reference):
    """
    Calculate the angle between v1 and v2.
    """
    angle = dot(v1, v2)/(norm(v1)*norm(v2))
    sign = sign(cross(v1, reference))
    return sign*angle
