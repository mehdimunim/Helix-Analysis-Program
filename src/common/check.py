import os
import sys

def check_argument(arguments):
    """
     Vérifie si le nom de fichier passé en argument existe.
    """
    if len(arguments) == 2:
        file_name = arguments[1]
    else:
        message = """
        ERROR: missing pdb filename as argument
        usage: %s file.pdb""" % (arguments[0])
        sys.exit(message)
    if not os.path.exists(file_name):
        sys.exit("ERROR: file %s does not seem to exist" % (file_name))

    return file_name
