from ..parser import parse
from .tpr import tpr
from ..dssp import dssp as dssp_mod
from ..axis import principal_axis


def trajectory_tpr(trajectory_file):
    """
    Get tpr list for each frame

    ---
    Parameters:
    trajectory_file: concatenation of pdb files

    ---
    Return:
    list_thetas: list of lists of tprs. 
    Each sub-list corresponds to a frame and each item of the sub-list to a helix

    """

    backbones = parse(trajectory_file)
    list_thetas = []
    for i, backbone in enumerate(backbones):
        dssp = dssp_mod.DSSP(backbone)
        list_helices = dssp.get_ca()
        thetas = []
        for helix in list_helices:
            orig, axis = principal_axis(helix)
            theta = tpr(helix, axis, orig)
            thetas.append(theta)
        list_thetas.append(thetas)
    return list_thetas