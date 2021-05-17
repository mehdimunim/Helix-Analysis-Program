import numpy


def principal_axis2(alpha_carbons):
    """
    Calculate principal inertia axis for the structure along with its geometrical center

    ---
    Parameters:
    alpha_carbons: alpha carbons of the structure

    ---
    Return:

    center: geometrical center of the structure
    axis_direction: direction of the axis

    """
    # alpha carbons coordinates as a numpy array
    coord = numpy.array(alpha_carbons, float)

    # get geometrical center
    center = numpy.mean(coord, 0)
    coord = coord - center

    # create inertia matrix and extract eigenvectors and values
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)

    # sort eigenvalues
    order = numpy.argsort(e_values)

    # axis1 is the principal axis with the greatest eigenvalue
    _, _, axis1 = e_vectors[:, order].transpose()

    axis_direction = axis1 / numpy.linalg.norm(axis1)


    return axis_direction
