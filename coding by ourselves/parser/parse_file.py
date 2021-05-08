

def getPdbList():
    """transforme le fichier pdb en une liste contenant la position de tous les
       atomes 'CA'
    """

    pdb = open("glut1.pdb", "r")

    # stocker x,y,z des résidues extraites
    res_coords = [] 

    """ charged_res = ["ARG", "HIS", "LYS", "ASP", "GLU"] """

    # Itérer sur les lignes dans pdb
    for line in pdb:
        #  Vérifier si la ligne commence par "ATOM"
        if line.startswith('ATOM'):
            # Itérer sur les lignes dans pdb
            atom_id = line[6:11].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_name = line[21:22].strip()
            residue_id = line[22:26].strip()
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
        
            for i in atom_id :
                if atom_name == 'CA' :
                    #res_coords.append([atom_id x, y, z])
                    res_coords.append([float(x), float(y), float(z)])
    return res_coords


def getPdbDict():
    """transforme le fichier pdb en un dictionnaire contenant une liste de tous les atomes
       pour chaque coordonnée de carbone 'CA'
    """

    pdb = open("glut1.pdb", "r")

    """ On stocke les résidus dans un dictionnaire de la manière suivante:
        les clés sont un tuple de 3 éléments (x,y,z) qui représentent les coordonnées du début du résidut
        les valeurs sont une liste de chaque atome dans ce résidu rangés par ordre d'apparition
    """
    #déclaration du dictionnaire qui comportera toutes les informations
    res_coords = {}

    #On déclare les variables qui serviront à stocker le début de chaque résidu:
    residueX, residueY, residueZ = 'none','none','none'

    # Itérer sur les lignes dans pdb
    for line in pdb:
        #  Vérifier si la ligne commence par "ATOM"
        if line.startswith('ATOM'):
            # Itérer sur les lignes dans pdb
            atom_id = line[6:11].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_name = line[21:22].strip()
            residue_id = line[22:26].strip()

            if atom_name == 'CA' :
                residueX = line[30:38].strip()
                residueY = line[38:46].strip()
                residueZ = line[46:54].strip()
                #on créé l'entrée définie par la clé (x,y,z) et on lui attribue une liste contenant le 1er atome
                #rencontré
                res_coords[(float(residueX),float(residueY),float(residueZ))] = ['CA']

            #sinon on ajoute l'atome à la liste du résidu en question
            #pour retrouver ce résidu, on le cherche dans le dictionnaire grâce à ses coordoonées
            #on fait attention à vérifier que cette liste existe bien
            elif residueX != 'none' and residueY !='none' and residueZ !='none':

                res_coords[(float(residueX),float(residueY),float(residueZ))].append(atom_name)
            else:
                residueX = line[30:38].strip()
                residueY = line[38:46].strip()
                residueZ = line[46:54].strip()
                res_coords[(float(residueX),float(residueY),float(residueZ))] = [atom_name]
    
    return res_coords





