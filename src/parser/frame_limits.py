def get_frames_limits(filename):
    """
    Get starts and ends of each frames a trajectory file
    Get start and end of the file for a static file
    ---
    Parameters:
    filename: trajectory or static pdb file

    ---
    Return:
    frames: list of tuple with (startline, endline) for each frame

    """
    frames = []
    with open(filename, "r") as pdbname:
        lines = pdbname.readlines()
        start = 0
        end = len(lines)
        for i, line in enumerate(lines):
            if line.startswith("MODEL"):
                start = i

            elif line.startswith("ENDMDL"):
                end = i
                frames.append((start, end))
    if len(frames) == 0:
        frames.append((start, end))
    return frames
