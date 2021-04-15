#
# Cleans the PDB file from irregularities
#

def check_irregularities(pdb):
    """
    Remove all proline residues from a pdb file
    Output: input pdb file without proline residues
    """
    with open(pdb + ".pdb", "r") as file_input:
        with open(pdb + "_cleaned.pdb", "w+") as output:
            for line in file_input:
                if line.startswith('ATOM'):
                    res_name = line[17:20].strip()
                    if res_name != 'PRO':
                        output.write(line)
    return pdb + "_cleaned.pdb"
