def check_irregularities():
    """
    Remove all proline residues from a pdb file
    Output: input pdb file without proline residues
    """
    with open("glut1.pdb", "r") as file_input:
        with open("glut1_removed.pdb", "w") as output:
            for line in file_input:
                if line.startswith('ATOM'):
                    res_name = line[17:20].strip()
                    if res_name != 'PRO':
                        output.write(line)
