

pdb = open("glut1.pdb", "r")

# stocker x,y,z des résiudes extraites
res_coords = [] 

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
    
        if atom_name == 'CA' :
            res_coords.append([x, y, z])

for item in res_coords:
        print (item)