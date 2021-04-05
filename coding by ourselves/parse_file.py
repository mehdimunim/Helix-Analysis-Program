

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
            res_coords[(residueX,residueY,residueZ)] = ['CA']

        #sinon on ajoute l'atome à la liste du résidu en question
        #pour retrouver ce résidu, on le cherche dans le dictionnaire grâce à ses coordoonées
        #on fait attention à vérifier que cette liste existe bien
        elif residueX != 'none' and residueY !='none' and residueZ !='none':

            res_coords[(residueX,residueY,residueZ)].append(atom_name)
        else:
            residueX = line[30:38].strip()
            residueY = line[38:46].strip()
            residueZ = line[46:54].strip()
            res_coords[(residueX,residueY,residueZ)] = [atom_name]

for item in res_coords:
        print (item," :",res_coords[item])