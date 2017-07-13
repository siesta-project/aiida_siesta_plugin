"""A setuptools based AiiDA plugin setup module.
Based on tutorial page:
http://aiida-core.readthedocs.io/en/latest/developer_guide/plugins/update_plugin.html

For packaging help:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

import json

from setuptools import find_packages, setup

if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(include_package_data=True, packages=find_packages(), **kwargs)
