#!/usr/bin/env runaiida
import os
import sys
import aiida_siesta.data.psf as psf
from aiida.cmdline.utils import decorators

#Check needed on the --stop-if-existing option

@decorators.with_dbenv()
def uploadfamily(*args):
    """
    Upload a new PSF-pseudopotential family.
    Returns the numbers of files found and the number of nodes uploaded.
    Call without parameters to get some help.
    """
    print(args)

    if not len(args) == 3 and not len(args) == 4:
        print(("Usage:"), file=sys.stderr)
        print(("runaiida uploadfamily.py FOLDER_NAME <group_name>  <group_description> "
                              "[OPTIONAL: --stop-if-existing]\n"), file=sys.stderr)
        sys.exit(1)

    folder = args[0]
    group_name = args[1]
    group_description = args[2]
    stop_if_existing = False

    if len(args) == 4:
        if args[3] == "--stop-if-existing":
            stop_if_existing = True
        else:
            print('Unknown directive: ' + args[3], file=sys.stderr)
            sys.exit(1)

    if (not os.path.isdir(folder)):
        print('Cannot find directory: ' + folder, file=sys.stderr)
        sys.exit(1)

    files_found, files_uploaded = psf.upload_psf_family(
        folder, group_name, group_description, stop_if_existing)

    print("PSF files found: {}. New files uploaded: {}".format(
        files_found, files_uploaded))


if __name__ == "__main__":
    #    uploadfamily(os.path.dirname(__file__), *sys.argv[1:])
    uploadfamily(*sys.argv[1:])
