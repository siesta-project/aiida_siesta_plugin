# -*- coding: utf-8 -*-
"""
Verdi command definition for PSML pseudopotentials data.
`verdi data psml`
"""
from __future__ import absolute_import
import click
import sys

from aiida.cmdline.commands.cmd_data import verdi_data
from aiida.cmdline.utils import decorators

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@verdi_data.group('psml', context_settings=CONTEXT_SETTINGS)
def psmldata():
    """PSML pseudos command line interface for SIESTA plugin."""
    pass


@psmldata.command()
@click.option('-n', '--name', prompt="Enter group name", help="Group name.")
@click.option(
    '-d',
    '--description',
    prompt="Enter group description",
    help="Brief description of the group.")
@click.option(
    '-S',
    '--stop-if-existing',
    is_flag=True,
    help=
    "OPTIONAL: Stop if some pseudo files are already registered in the database."
)
@click.argument(
    'directory',
    type=click.Path(exists=True, file_okay=False, resolve_path=True))
@decorators.with_dbenv()
def uploadfamily(name, description, stop_if_existing, directory):
    """
    Upload a collection of *.psml files from specified directory
    to a new pseudopotential family.

    Returns the numbers of files found and the number of nodes uploaded.
    """

    import aiida_siesta.data.psml as psml

    files_found, files_uploaded = psml.upload_psml_family(
        directory, name, description, stop_if_existing)

    click.echo("PSML files found: {}. New files uploaded: {}".format(
        files_found, files_uploaded))


@psmldata.command()
@click.option(
    '-e',
    '--element',
    multiple=True,
    help="OPTIONAL: Filter the families only to those containing "
    "a pseudo for each of the specified elements.")
@click.option(
    '-D',
    '--with-description',
    is_flag=True,
    help="OPTIONAL: Print families\' description.")
@decorators.with_dbenv()
def listfamilies(element, with_description):
    """
    Print on screen the list of installed PSML-pseudo families.
    """
    from aiida.plugins import DataFactory
    from aiida_siesta.data.psml import PSMLGROUP_TYPE

    PsmlData = DataFactory('siesta.psml')
    from aiida.orm import QueryBuilder
    from aiida.orm import Group
    qb = QueryBuilder()
    qb.append(PsmlData, tag='psmldata')

    if element:
        qb.add_filter(PsmlData, {'attributes.element': {'in': element}})

    qb.append(
        Group,
        group_of='psmldata',
        tag='group',
        project=["label", "description"],
        filters={"type_string": {
            '==': PSMLGROUP_TYPE
        }})

    qb.distinct()
    if qb.count() > 0:
        for res in qb.dict():
            group_name = res.get("group").get("label")
            group_desc = res.get("group").get("description")
            qb = QueryBuilder()
            qb.append(
                Group,
                tag='thisgroup',
                filters={"label": {
                    'like': group_name
                }})
            qb.append(PsmlData, project=["id"], member_of='thisgroup')

            if with_description:
                description_string = ": {}".format(group_desc)
            else:
                description_string = ""

            click.echo("* {} [{} pseudos]{}".format(group_name, qb.count(),
                                                    description_string))

    else:
        click.echo("No valid PSML pseudopotential family found.", err=True)


@psmldata.command()
@click.argument('family')
@click.argument(
    'directory',
    type=click.Path(exists=False, file_okay=False, resolve_path=True))
@decorators.with_dbenv()
def exportfamily(family, directory):
    """
    Export a pseudopotential family into a new directory.
    """
    import os
    import io
    from aiida.common.exceptions import NotExistent
    from aiida.plugins import DataFactory

    PsmlData = DataFactory('siesta.psml')
    try:
        group = PsmlData.get_psml_group(family)
    except NotExistent:
        click.echo("PSML family {} not found".format(family), err=True)

    try:
        os.makedirs(directory)
        for pseudo in group.nodes:
            dest_path = os.path.join(directory, pseudo.filename)
            # with open(dest_path, 'w') as dest:
            with io.open(dest_path, 'w', encoding='utf8') as handle:
                # with pseudo._get_folder_pathsubfolder.open(
                #         pseudo.filename) as source:
                #     dest.write(source.read())
                handle.write(pseudo.get_content())
    except OSError:
        click.echo(
            "Destination directory {} exists; aborted".format(directory),
            err=True)
