#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os
import sys
from aiida import load_dbenv
load_dbenv()
import aiida_siesta.data.psf as psf

def uploadfamily(*args):
    """
    Upload a new PSF-pseudopotential family.
    Returns the numbers of files found and the number of nodes uploaded.
    Call without parameters to get some help.
    """
    print args
    
    if not len(args) == 3 and not len(args) == 4:
        print >> sys.stderr, ("Usage:")
        print >> sys.stderr, ("./uploadfamily.py FOLDER_NAME <group_name>  <group_description> "
                              "[OPTIONAL: --stop-if-existing]\n")
        sys.exit(1)

    folder = args[0]
    group_name = args[1]
    group_description = args[2]
    stop_if_existing = False

    if len(args) == 4:
        if args[3] == "--stop-if-existing":
            stop_if_existing = True
        else:
            print >> sys.stderr, 'Unknown directive: ' + args[3]
            sys.exit(1)

    if (not os.path.isdir(folder)):
        print >> sys.stderr, 'Cannot find directory: ' + folder
        sys.exit(1)

    files_found, files_uploaded = psf.upload_psf_family(
        folder, group_name, group_description, stop_if_existing)

    print "PSF files found: {}. New files uploaded: {}".format(
        files_found, files_uploaded)


if __name__ == "__main__":
    #    uploadfamily(os.path.dirname(__file__), *sys.argv[1:])
    uploadfamily(*sys.argv[1:])
