"""
28-April-2020
This module manages the LUA in the local repository. Its Similar to Pseudoptenial psf/psml.
"""

import io
from aiida.common.utils import classproperty
from aiida.common.files import md5_file
from aiida.orm.nodes import SinglefileData


def parse_lua(fname, check_filename=True):
    """
    Try to get relevant information from the LUA. For the moment, only the
    element name.
    If check_filename is True, raise a ParsingError exception if the filename
    does not start with the element name.
    """
    import os

    from aiida.common.exceptions import ParsingError
    # from aiida.common import AIIDA_LOGGER
    from aiida.orm.nodes.data.structure import _valid_symbols

    parsed_data = {}

    try:
        lua_contents = fname.read().split()
        fname = fname.name
    except AttributeError:
        with io.open(fname, encoding='utf8') as fil:
            lua_contents = fil.read().split()

    # Parse the element
    #element = None
    element = fname
    for element in lua_contents:
        break

    # Only first letter capitalized!
    if element is None:
        raise ParsingError("Unable to find the element of lua {}".format(fname))
    element = element.capitalize()
#    if element not in _valid_symbols:
#        raise ParsingError("Unknown element symbol {} for file {}".format(element, fname))

#    if check_filename:
#        if not os.path.basename(fname).lower().startswith(element.lower()):
#            raise ParsingError(
#                "Filename {0} was recognized for element "
#                "{1}, but the filename does not start "
#                "with {1}".format(fname, element)
#            )

    parsed_data['element'] = element

    return parsed_data


class LUAData(SinglefileData):
    """
    Function not yet documented.
    """

    @classmethod
    def get_or_create(cls, filename, use_first=False, store_lua=True):
        """
        Pass the same parameter of the init; if a file with the same md5
        is found, that LUAData is returned.

        :param filename: an absolute filename on disk
        :param use_first: if False (default), raise an exception if more than \
                one potential is found.\
                If it is True, instead, use the first available lua.
        :param bool store_lua: If false, the LUAData objects are not stored in
                the database. default=True.
        :return (lua, created): where lua is the LUAData object, and create is either\
            True if the object was created, or False if the object was retrieved\
            from the DB.
        """
        import os

        if not os.path.abspath(filename):
            raise ValueError("filename must be an absolute path")

        md5 = md5_file(filename)

        luafile = cls.from_md5(md5)
        if not luafile:
            instance = cls(file=filename)
            if store_lua:
                instance.store()
            return (instance, True)

        if len(luafile) > 1:
            if use_first:
                return (luafile[0], False)

            raise ValueError(
                "More than one copy of lua file "
                "with the same MD5 has been found in the "
                "DB. pks={}".format(",".join([str(i.pk) for i in luafile]))
            )

        return (luafile[0], False)



    def store(self, *args, **kwargs):  # pylint: disable=arguments-differ
        """
        Store the node, reparsing the file so that the md5 and the element
        are correctly reset.
        """
        from aiida.common.exceptions import ParsingError
        from aiida.common.files import md5_from_filelike

        if self.is_stored:
            return self

        with self.open(mode='r') as handle:
            parsed_data = parse_lua(handle)

        # Open in binary mode which is required for generating the md5 checksum
        with self.open(mode='rb') as handle:
            md5sum = md5_from_filelike(handle)

        try:
            element = parsed_data['element']
        except KeyError:
            raise ParsingError("No 'element' parsed in the LUA file {};" " unable to store".format(self.filename))

        self.set_attribute('element', str(element))
        self.set_attribute('md5', md5sum)

        return super(LUAData, self).store(*args, **kwargs)

    @classmethod
    def from_md5(cls, md5):
        """
        Return a list of all LUA files that match a given MD5 hash.

        Note that the hash has to be stored in a _md5 attribute, otherwise
        the luafile will not be found.
        """
        from aiida.orm.querybuilder import QueryBuilder
        qb = QueryBuilder()
        qb.append(cls, filters={'attributes.md5': {'==': md5}})
        return [_ for [_] in qb.all()]

    def set_file(self, filename):  # pylint: disable=arguments-differ
        """
        I pre-parse the file to store the attributes.
        """
        from aiida.common.exceptions import ParsingError

        parsed_data = parse_lua(filename)
        md5sum = md5_file(filename)

        try:
            element = parsed_data['element']
        except KeyError:
            raise ParsingError("No 'element' parsed in the LUA file {};" " unable to store".format(self.filename))

        super(LUAData, self).set_file(filename)

        # self._set_attr('element', str(element))
        # self._set_attr('md5', md5sum)
        self.set_attribute('element', str(element))
        self.set_attribute('md5', md5sum)


    @property
    def element(self):
        return self.get_attribute('element', None)

    @property
    def md5sum(self):
        return self.get_attribute('md5', None)

    def _validate(self):
        from aiida.common.exceptions import ValidationError
        from aiida.common.files import md5_from_filelike

        super(LUAData, self)._validate()

        with self.open(mode='r') as handle:
            parsed_data = parse_lua(handle)

        # Open in binary mode which is required for generating the md5 checksum
        with self.open(mode='rb') as handle:
            md5 = md5_from_filelike(handle)

        # This is erroneous exception,
        # as it is in the `upf` module oin `aiida_core`
        try:
            element = parsed_data['element']
        except KeyError:
            raise ValidationError("No 'element' could be parsed in the LUA " "file {}".format(self.filename))

        try:
            attr_element = self.get_attribute('element')
        except AttributeError:
            raise ValidationError("attribute 'element' not set.")

        try:
            attr_md5 = self.get_attribute('md5')
        except AttributeError:
            raise ValidationError("attribute 'md5' not set.")

        if attr_element != element:
            raise ValidationError(
                "Attribute 'element' says '{}' but '{}' was "
                "parsed instead.".format(attr_element, element)
            )

        if attr_md5 != md5:
            raise ValidationError("Attribute 'md5' says '{}' but '{}' was " "parsed instead.".format(attr_md5, md5))



