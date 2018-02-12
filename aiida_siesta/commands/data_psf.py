# -*- coding: utf-8 -*-
"""
Verdi command definition for PSF pseudopotentials data.
`verdi data psf`
"""
import click
import sys

from aiida.cmdline.commands import data_cmd
from aiida.backends.utils import is_dbenv_loaded, load_dbenv

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@data_cmd.group('psf', context_settings=CONTEXT_SETTINGS)
def psfdata():
    """PSF pseudos command line interface for SIESTA plugin."""
    pass


@psfdata.command()
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
def uploadfamily(name, description, stop_if_existing, directory):
    """
    Upload a collection of *.psf files from specified directory
    to a new pseudopotential family.

    Returns the numbers of files found and the number of nodes uploaded.
    """
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

    import aiida_siesta.data.psf as psf

    files_found, files_uploaded = psf.upload_psf_family(
        directory, name, description, stop_if_existing)

    click.echo("PSF files found: {}. New files uploaded: {}".format(
        files_found, files_uploaded))


@psfdata.command()
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
def listfamilies(element, with_description):
    """
    Print on screen the list of installed PSF-pseudo families.
    """
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

    from aiida.orm import DataFactory
    from aiida_siesta.data.psf import PSFGROUP_TYPE

    PsfData = DataFactory('siesta.psf')
    from aiida.orm.querybuilder import QueryBuilder
    from aiida.orm.group import Group
    qb = QueryBuilder()
    qb.append(PsfData, tag='psfdata')

    if element:
        qb.add_filter(PsfData, {'attributes.element': {'in': element}})

    qb.append(
        Group,
        group_of='psfdata',
        tag='group',
        project=["name", "description"],
        filters={
            "type": {
                '==': PSFGROUP_TYPE
            }
        })

    qb.distinct()
    if qb.count() > 0:
        for res in qb.dict():
            group_name = res.get("group").get("name")
            group_desc = res.get("group").get("description")
            qb = QueryBuilder()
            qb.append(
                Group, tag='thisgroup', filters={
                    "name": {
                        'like': group_name
                    }
                })
            qb.append(PsfData, project=["id"], member_of='thisgroup')

            if with_description:
                description_string = ": {}".format(group_desc)
            else:
                description_string = ""

            click.echo("* {} [{} pseudos]{}".format(group_name,
                                                    qb.count(),
                                                    description_string))

    else:
        click.echo("No valid PSF pseudopotential family found.", err=True)


@psfdata.command()
@click.argument('family')
@click.argument(
    'directory',
    type=click.Path(exists=False, file_okay=False, resolve_path=True))
def exportfamily(family, directory):
    """
    Export a pseudopotential family into a new directory.
    """
    from aiida import is_dbenv_loaded, load_dbenv
    if not is_dbenv_loaded():
        load_dbenv()

    import os
    from aiida.common.exceptions import NotExistent
    from aiida.orm import DataFactory

    PsfData = DataFactory('siesta.psf')
    try:
        group = PsfData.get_psf_group(family)
    except NotExistent:
        click.echo("PSF family {} not found".format(family), err=True)

    try:
        os.makedirs(directory)
        for pseudo in group.nodes:
            dest_path = os.path.join(directory, pseudo.filename)
            with open(dest_path, 'w') as dest:
                with pseudo._get_folder_pathsubfolder.open(
                        pseudo.filename) as source:
                    dest.write(source.read())
    except OSError:
        click.echo(
            "Destination directory {} exists; aborted".format(directory),
            err=True)
