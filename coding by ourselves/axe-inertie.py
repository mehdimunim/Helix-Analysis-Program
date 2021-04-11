
"""
Computes principal axes from a PDB file.
OL.

tape python3 axe_inertie.py  glut1.pdb  sur le terminal
"""

import sys
import os.path
import numpy




def read_pdb_xyz(pdb_name):
    """
    Reads atomic coordinates of C-alpha atoms from a .pdb file.
    Parameters
    ----------
    pdb_name : str
        Name of pdb file.
    Returns
    -------
    array of atomic coordinates
        [[x1 y1 z1]
         [x2 y2 z2]
         [.. .. ..]
         [xn yn zn]]
    """
    xyz = []
    with open(pdb_name, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                # extract x, y, z coordinates for carbon alpha atoms
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                if line[12:16].strip() == "CA":
                    xyz.append([x, y, z])
    return xyz


def check_argument(arguments):
    """
    Check if filename passed as argument exists. 
    Parameters
    """
   
    if len(arguments) == 2:
        file_name = arguments[1]
    else:
        message = """
        ERROR: missing pdb filename as argument
        usage: %s file.pdb""" %(arguments[0])
        sys.exit(message)

    # check if argument is an existing file
    if not os.path.exists(file_name):
        sys.exit("ERROR: file %s does not seem to exist" %(file_name))

    return file_name



if __name__ == '__main__':

    # check if argument file exists
    pdb_name = check_argument(sys.argv)


    xyz = read_pdb_xyz(pdb_name)
    print("%d CA atomes found in %s" %(len(xyz), pdb_name))

    #create coordinates array
    coord = numpy.array(xyz, float)

    # compute geometric center
    center = numpy.mean(coord, 0)
    print("Coordinates of the geometric center:\n", center)

    # center with geometric center
    coord = coord - center

    # compute principal axis matrix
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)
 
    order = numpy.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
  
    print("\nFirst principal axis")
    print("coordinates: ", axis1)
   

   