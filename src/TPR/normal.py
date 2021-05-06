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
    # vector from axis normal_origin to CA
    vector = point - axis_origin

    # projection of vector on axis
    projection = np.dot(vector, axis)*axis_direction

    # the normal vector from normal_origin to point
    normal_direction = projection - vector

    # normalized normal_direction
    normal_direction = normal_direction/np.linalg.norm(normal_direction)

    # opposite of the projection of vector on normal 
    projection2 = - np.dot(normal_direction, vector)*normal_direction

    normal_origin = point + projection2

    return normal_direction, normal_origin
