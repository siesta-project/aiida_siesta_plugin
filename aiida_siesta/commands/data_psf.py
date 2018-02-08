# -*- coding: utf-8 -*-
"""
Verdi command definition for PSF pseudopotentials data.
`verdi data psf`
"""
import click
import sys

from aiida.cmdline.commands import data_cmd
from aiida.backends.utils import is_dbenv_loaded, load_dbenv

@data_cmd.group('psf')
def psfdata():
    """Command line interface for template plugin"""
    pass


@psfdata.command()
# def listfamilies(*args):
def listfamilies():
    """
    Print on screen the list of psf families installed
    """
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

    # TODO Adapt arguments parsing for new Click interface
    # note that the following command requires that the psfdata has a
    # key called element. As such, it is not well separated.
    # import argparse

    # parser = argparse.ArgumentParser(
    #     prog=self.get_full_command_name(),
    #     description='List AiiDA psf families.')
    # parser.add_argument(
    #     '-e',
    #     '--element',
    #     nargs='+',
    #     type=str,
    #     default=None,
    #     help="Filter the families only to those containing "
    #     "a pseudo for each of the specified elements")
    # parser.add_argument(
    #     '-d',
    #     '--with-description',
    #     dest='with_description',
    #     action='store_true',
    #     help="Show also the description for the PSF family")
    # parser.set_defaults(with_description=False)

    # args = list(args)
    # parsed_args = parser.parse_args(args)

    from aiida.orm import DataFactory
    from aiida_siesta.data.psf import PSFGROUP_TYPE

    PsfData = DataFactory('siesta.psf')
    from aiida.orm.querybuilder import QueryBuilder
    from aiida.orm.group import Group
    qb = QueryBuilder()
    qb.append(PsfData, tag='psfdata')
    # if parsed_args.element is not None:
    #     qb.add_filter(PsfData,
    #                     {'attributes.element': {
    #                         'in': parsed_args.element
    #                     }})
    qb.append(
        Group,
        group_of='psfdata',
        tag='group',
        project=["name", "description"],
        filters={"type": {
            '==': PSFGROUP_TYPE
        }})

    qb.distinct()
    if qb.count() > 0:
        for res in qb.dict():
            group_name = res.get("group").get("name")
            group_desc = res.get("group").get("description")
            qb = QueryBuilder()
            qb.append(
                Group,
                tag='thisgroup',
                filters={"name": {
                    'like': group_name
                }})
            qb.append(PsfData, project=["id"], member_of='thisgroup')

            # if parsed_args.with_description:
            #     description_string = ": {}".format(group_desc)
            # else:
            #     description_string = ""
            description_string = ""

            print "* {} [{} pseudos]{}".format(group_name,
                                                qb.count(),
                                                description_string)

    else:
        print "No valid PSF pseudopotential family found."

##################
# from aiida.cmdline.baseclass import VerdiCommandWithSubcommands
# from aiida.cmdline.commands.data import Importable


# class _Psf(VerdiCommandWithSubcommands, Importable):
#     """
#     Setup and manage psf to be used

#     This command allows to list and configure psf pseudos.
#     """

#     def __init__(self):
#         """
#         A dictionary with valid commands and functions to be called.
#         """
#         if not is_dbenv_loaded():
#             load_dbenv()
#         # from aiida.orm.data.upf import UpfData
#         from aiida_siesta.data.psf import PsfData

#         self.dataclass = PsfData
#         self.valid_subcommands = {
#             'uploadfamily': (self.uploadfamily, self.complete_auto),
#             'listfamilies': (self.listfamilies, self.complete_none),
#             'import': (self.importfile, self.complete_none),
#             'exportfamily': (self.exportfamily, self.complete_auto)
#         }

#     def uploadfamily(self, *args):
#         """
#         Upload a new pseudopotential family.

#         Returns the numbers of files found and the number of nodes uploaded.

#         Call without parameters to get some help.
#         """
#         import os.path

#         if not len(args) == 3 and not len(args) == 4:
#             print >> sys.stderr, (
#                 "After 'psf uploadfamily' there should be three "
#                 "arguments:")
#             print >> sys.stderr, ("folder, group_name, group_description "
#                                   "[OPTIONAL: --stop-if-existing]\n")
#             sys.exit(1)

#         folder = os.path.abspath(args[0])
#         group_name = args[1]
#         group_description = args[2]
#         stop_if_existing = False

#         if len(args) == 4:
#             if args[3] == "--stop-if-existing":
#                 stop_if_existing = True
#             else:
#                 print >> sys.stderr, 'Unknown directive: ' + args[3]
#                 sys.exit(1)

#         if (not os.path.isdir(folder)):
#             print >> sys.stderr, 'Cannot find directory: ' + folder
#             sys.exit(1)

#         # import aiida.orm.data.upf as upf
#         import aiida_siesta.data.psf as psf

#         files_found, files_uploaded = psf.upload_psf_family(
#             folder, group_name, group_description, stop_if_existing)

#         print "PSF files found: {}. New files uploaded: {}".format(
#             files_found, files_uploaded)

#     def listfamilies(self, *args):
#         """
#         Print on screen the list of psf families installed
#         """
#         # note that the following command requires that the psfdata has a
#         # key called element. As such, it is not well separated.
#         import argparse

#         parser = argparse.ArgumentParser(
#             prog=self.get_full_command_name(),
#             description='List AiiDA psf families.')
#         parser.add_argument(
#             '-e',
#             '--element',
#             nargs='+',
#             type=str,
#             default=None,
#             help="Filter the families only to those containing "
#             "a pseudo for each of the specified elements")
#         parser.add_argument(
#             '-d',
#             '--with-description',
#             dest='with_description',
#             action='store_true',
#             help="Show also the description for the PSF family")
#         parser.set_defaults(with_description=False)

#         args = list(args)
#         parsed_args = parser.parse_args(args)

#         from aiida.orm import DataFactory
#         # from aiida.orm.data.upf import UPFGROUP_TYPE
#         from aiida_siesta.data.psf import PSFGROUP_TYPE

#         PsfData = DataFactory('siesta.psf')
#         from aiida.orm.querybuilder import QueryBuilder
#         from aiida.orm.group import Group
#         qb = QueryBuilder()
#         qb.append(PsfData, tag='psfdata')
#         if parsed_args.element is not None:
#             qb.add_filter(PsfData,
#                           {'attributes.element': {
#                               'in': parsed_args.element
#                           }})
#         qb.append(
#             Group,
#             group_of='psfdata',
#             tag='group',
#             project=["name", "description"],
#             filters={"type": {
#                 '==': PSFGROUP_TYPE
#             }})

#         qb.distinct()
#         if qb.count() > 0:
#             for res in qb.dict():
#                 group_name = res.get("group").get("name")
#                 group_desc = res.get("group").get("description")
#                 qb = QueryBuilder()
#                 qb.append(
#                     Group,
#                     tag='thisgroup',
#                     filters={"name": {
#                         'like': group_name
#                     }})
#                 qb.append(PsfData, project=["id"], member_of='thisgroup')

#                 if parsed_args.with_description:
#                     description_string = ": {}".format(group_desc)
#                 else:
#                     description_string = ""

#                 print "* {} [{} pseudos]{}".format(group_name,
#                                                    qb.count(),
#                                                    description_string)

#         else:
#             print "No valid PSF pseudopotential family found."

#     def exportfamily(self, *args):
#         """
#         Export a pseudopotential family into a folder.
#         Call without parameters to get some help.
#         """
#         import os
#         from aiida.common.exceptions import NotExistent
#         from aiida.orm import DataFactory

#         if not len(args) == 2:
#             print >> sys.stderr, ("After 'psf export' there should be two "
#                                   "arguments:")
#             print >> sys.stderr, ("folder, psf_family_name\n")
#             sys.exit(1)

#         folder = os.path.abspath(args[0])
#         group_name = args[1]

#         PsfData = DataFactory('siesta.psf')
#         try:
#             group = PsfData.get_psf_group(group_name)
#         except NotExistent:
#             print >> sys.stderr, ("psf family {} not found".format(group_name))

#         for u in group.nodes:
#             dest_path = os.path.join(folder, u.filename)
#             if not os.path.isfile(dest_path):
#                 with open(dest_path, 'w') as dest:
#                     with u._get_folder_pathsubfolder.open(
#                             u.filename) as source:
#                         dest.write(source.read())
#             else:
#                 print >> sys.stdout, ("File {} is already present in the "
#                                       "destination folder".format(u.filename))

#     def _import_psf(self, filename, **kwargs):
#         """
#         Importer from PSF.
#         """
#         try:
#             node, _ = self.dataclass.get_or_create(filename)
#             print node
#         except ValueError as e:
#             print e
