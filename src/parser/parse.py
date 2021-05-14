#
# Gets the backbone atoms
#

from .frame_limits import get_frames_limits
from .parse_structure import parse_structure


def parse(filename, getTraj=False):
    """
    Parse a pdb file to get all backbone atoms

    ---
    Parameters:
    filename: trajectory or static pdb file

    ---
    Returns:
    if trajectory, list of backbones for each frames
    else, backbone of the molecule

    """
    limits = get_frames_limits(filename)
    backbones = []
    for limit in limits:
        backbone = parse_structure(filename, limit[0], limit[1])
        backbones.append(backbone)
    if getTraj:
        if (len(backbones) == 1):
            return (backbones[0], False)
        return (backbones, True)
    else:
        if (len(backbones) == 1):
            return backbones[0]
        return backbones
