# aiida_siesta_plugin

This repository contains the files that implement the AiiDA Siesta
plugin (input, parser, and tests) as a stand-alone package.

## DEVELOPMENT AND USAGE guidelines:

# Local development

Use your own installation of the database, and your own computers and codes.

You need to work on a specially patched copy of aiida_core that works
around a few issues relevant to the Siesta plugin:

  git clone https://github.com/vdikan/aiida_core
  cd aiida_core
  git checkout v0.9_siesta_patch

Install aiida with the '-e' option (and optionally install the
dependencies for extra utilities such as pymatgen, cif, etc, and the
generation of docs)

	   cd aiida_core
	   pip install -e .
	   pip install -e .[docs,atomic_tools]

The first time you will need to configure AiiDA as explained in the manual.

Install the plugin by executing, from the top level of the plugin directory:

	pip install -e .

# Docker-based usage and automated testing

The Docker framework is being updated to the new plugin architecture.

