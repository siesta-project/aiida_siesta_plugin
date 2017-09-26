"""A setuptools based AiiDA plugin setup module.
Based on tutorial page:
http://aiida-core.readthedocs.io/en/latest/developer_guide/plugins/update_plugin.html

For packaging help:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

import codecs
import json
import os

from setuptools import find_packages, setup

THIS_LOC = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(THIS_LOC, *parts), "rb", "utf-8") as f:
        return f.read()


if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(
        include_package_data=True,
        packages=find_packages(),
        long_description=read('README.rst'),
        reentry_register=True,
        **kwargs)
