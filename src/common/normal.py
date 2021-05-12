import numpy as np


def normal(axis_direction, axis_origin, point):
    """
    Calculate the normal to axis intercepting point
    ---
    Parameters:
    axis_direction: normal_direction of the axis
    axis_origin: starting point of the axis
    point: intercepts the normal

    ---
    Return:
    normal_direction: steering vector of normal 
    normal_origin: starting point of normal

    """
    # transform input into numpy arrays
    axis_direction = np.array(axis_direction, float)
    axis_origin = np.array(axis_origin, float)
    point = np.array(point, float)

    # vector from axis normal_origin to point
    vector = point - axis_origin

    # projection of vector on axis
    projection = np.dot(vector, axis_direction)*axis_direction

    # the normal vector from normal_origin to point
    normal_direction = vector - projection

    # normalized normal_direction
    normal_direction = normal_direction/np.linalg.norm(normal_direction)

    # opposite of the projection of vector on normal
    projection2 = - np.dot(normal_direction, vector)*normal_direction

    normal_origin = point + projection2

    return normal_direction, normal_origin


def test_normal():
    import math

    # Test when point is after axis_origin

    a = math.pi/6
    axis_origin = [0, 0, 0]
    axis_direction = [1, 0, 0]
    CA = [math.cos(a), math.sin(a), 0]

    t = 10**(-5)

    dir, orig = normal(axis_direction, axis_origin, CA)
    # dir should be [0, 1, 0]
    print(dir)
    assert abs(dir[1]-1) < t, "wrong direction"
    # orig should be [cos(a), 0, 0]
    print(orig)
    assert abs(orig[0] - math.cos(a) < t), "wrong origin"

    # Test when point is before axis_origin

    axis_origin = [3, 0, 0]
    # dir should be [0, 1, 0]
    print(dir)
    assert abs(dir[1]-1) < t, "wrong direction"

    # orig should be [cos(a), 0, 0]
    print(orig)
    assert abs(orig[0] - math.cos(a) < t), "wrong origin"


# test_normal()
