"""
Ce fichier calcule les axes d'inerties principaux d'une proteine à partir d'un fichier PDB.
Ici, les deux axes possédant les deux grandes valeurs propres. 
Calcule l'angle entre ces deux axes. 

Et afficher le graphe de variation de l'angle tout au long de la trajectiore
    --> ATTENTION : il faut prendre un fichier possédant une trajectoire : TSPO_traj.pdb 

Pour faire tourner le programme : taper "python3 trajectory.py TSPO_traj.pdb" dans le terminal 

AVANT : intaller module cogent3
        installer MDAnalysis
        installer biopython=1.68 

        
A ADAPTER aux hélices de la proteine 
"""

import sys
import os.path
import numpy
import math
import MDAnalysis as mda
import matplotlib.pyplot as plt
from ..parser import parse
from ..common import angle

from cogent3.maths.stats.special import fix_rounding_error


def inertia_helix(helix):
    """
    Calculate angles between inertia axes 1 and 2 of helix

    ---
    Parameters:
    helix : list of a-carbons of the helix

    ---
    Return:
    the angle between inertia axes in °

    """

    Rgyr = []
    # Création d'un tableau de coordonnées
    coord = numpy.array(helix, float)

    # Calcul les coordonnées du centre géométrique
    center = numpy.mean(coord, 0)
    coord = coord - center

    # Création de la matrice d'inertie et extraction des valeurs et vecteurs propres
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)

    # Ordonner les valeurs propres
    order = numpy.argsort(e_values)

    # axis1 est l'axe principal avec la plus grande valeur propre (val1)
    # axis2 est l'axe principal avec la deuxième plus grande valeur propre (val2).
    _, axis2, axis1 = e_vectors[:, order].transpose()

    # Calcul l'angle entre les l'axe1 et l'axe2
    """calcule l'angle entre deux vecteurs v1 et v2 qui sont des objets numpy.array.
    renvoie un flottant contenant l'angle en radians à convertir en degrée.
    """
    angle = angle(axis1, axis2)
    return math.degrees(angle)


def inertia(list_helices, pdb_name):
    """
    Calculate angles between inertia axes 1 and 2 of helix

    ---
    Output:
    Save it to output/inertia_MOLECULE

    """
    molecule_name = pdb_name.split("/")[1][:-4]

    angles = []
    for helix in list_helices:
        angle = inertia_helix(helix)
        angles.append(angle)

    angles = numpy.array(angles)
    ind = [i for i in range(len(list_helices))]

    plt.xlabel("#Helix")
    plt.ylabel("Angle (°)")
    plt.title("Angle between the first two inertia axes")
    plt.scatter(ind, angles)
    plt.savefig("output/inertia_" + molecule_name + ".png")


def inertia_traj(pdb_name):
    """
    Calculate angles between inertia axes 1 and 2

    ---
    Ouput:
    Save file in "ouput/inertia_traj_MOLECULE.png"

    """
    molecule_name = pdb_name.split("/")[1][:-4]
    backbone = parse(pdb_name)

    xyz = backbone["alpha_carbon"]
    u = mda.Universe(pdb_name)

    Rgyr = []
    for ts in u.trajectory:
        print(ts)
        # Création d'un tableau de coordonnées
        coord = numpy.array(xyz, float)

        # Calcul les coordonnées du centre géométrique
        center = numpy.mean(coord, 0)
        coord = coord - center

        # Création de la matrice d'inertie et extraction des valeurs et vecteurs propres
        inertia = numpy.dot(coord.transpose(), coord)
        e_values, e_vectors = numpy.linalg.eig(inertia)

        # Ordonner les valeurs propres
        order = numpy.argsort(e_values)

        # axis1 est l'axe principal avec la plus grande valeur propre (val1)
        # axis2 est l'axe principal avec la deuxième plus grande valeur propre (val2).
        _, axis2, axis1 = e_vectors[:, order].transpose()

        # Calcul l'angle entre les l'axe1 et l'axe2
        """calcule l'angle entre deux vecteurs v1 et v2 qui sont des objets numpy.array.
        renvoie un flottant contenant l'angle en radians à convertir en degrée.
        """
        angle = angle(axis1, axis2)

        Rgyr.append((u.trajectory.time, math.degrees(angle)))

    Rgyr = numpy.array(Rgyr)
    ax = plt.subplot(111)
    ax.plot(Rgyr[:, 0], Rgyr[:, 1], 'r--', lw=2, label=r"$R_G$")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"Angle between the first two inertia axes (°)")
    ax.figure.savefig("output/inertia_traj_" + molecule_name + ".png")
