"""
Ce fichier calcule les axes d'inerties principaux d'une proteine à partir d'un fichier PDB.
Ici, les deux axes possédant les deux grandes valeurs propres. 
Calcule l'angle entre ces deux axes. 


Pour faire tourner le programme : taper "python3 axe-inertie.py glut1.pdb" dansn le terminal 

AVANT : intaller module cogent3

A ADAPTER aux hélices de la proteine 
"""

import sys
import os.path
import numpy
import math

from cogent3.maths.stats.special import fix_rounding_error



def read_pdb_xyz(pdb_name):
    """
    Fonction qui retourne les coordonnées atomiques des atomes C-alpha à partir d'un fichier .pdb 
    sou forme de tableau de coordonnées
        [[x1 y1 z1]
         [x2 y2 z2]
         [.. .. ..]
         [xn yn zn]]
    """
    xyz = []
    with open(pdb_name, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                if line[12:16].strip() == "CA":
                    xyz.append([x, y, z])
    return xyz


def check_argument(arguments):
    """
     Vérifie si le nom de fichier passé en argument existe.
    """
    if len(arguments) == 2:
        file_name = arguments[1]
    else:
        message = """
        ERROR: missing pdb filename as argument
        usage: %s file.pdb""" %(arguments[0])
        sys.exit(message)
    if not os.path.exists(file_name):
        sys.exit("ERROR: file %s does not seem to exist" %(file_name))

    return file_name


def scalar(v1,v2):
    return sum(v1*v2)


if __name__ == '__main__':
    
    pdb_name = check_argument(sys.argv)

    xyz = read_pdb_xyz(pdb_name)

    #Création d'un tableau de coordonnées
    coord = numpy.array(xyz, float)

    # Calcul les coordonnées du centre géométrique
    center = numpy.mean(coord, 0)
    coord = coord - center

    # Création de la matrice d'inertie et extraction des valeurs et vecteurs propres
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)
 
    # Ordonner les valeurs propres 
    order = numpy.argsort(e_values)

    #axes1 est l'axe principal avec la plus grande valeur propre (val1)
    #axes2 est l'axe principal avec la deuxième plus grande valeur propre (val2).
    val3, val2, val1 = e_values[order]
    axes3, axes2, axes1 = e_vectors[:, order].transpose()


    #Calcul l'angle entre les l'axe1 et l'axe2
    """calcule l'angle entre deux vecteurs v1 et v2 qui sont des objets numpy.array.
    renvoie un flottant contenant l'angle en radians à convertir en degrée.
    """
    unit_vector_1 = axes1 / numpy.linalg.norm(axes1)
    unit_vector_2 = axes2 / numpy.linalg.norm(axes2)
    dot_product = numpy.dot(unit_vector_1, unit_vector_2)
    angle = numpy.arccos(dot_product)


    #Affichage :

    print("\n Axe d'inertie 1 : ")
    print("Coordonnées: ", axes1)
    print("Valeur propre :", val1)

    print("\n Axe d'inertie 2 :")
    print("Coordonnées:", axes2)
    print("Valeur propre :", val2)


    print("\n Angle entre les deux axes (en degré): " , math.degrees(angle) )

