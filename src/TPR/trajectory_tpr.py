
def trajectory_tpr(trajectory_file):
    """
    Get tpr list for each frame

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    tpr: list of list of tprs

    """

    for frame in list_frames
        list_helices = dssp.get_ca()

        for helix in list_helices:
            orig, axis = principal_axis(helix)
            theta = tpr(helix, axis, orig)
